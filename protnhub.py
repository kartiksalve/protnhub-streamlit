import streamlit as st
import requests
import networkx as nx
import plotly.graph_objects as go
import kaleido  # Ensure kaleido is installed (pip install kaleido)
import openai

# ------------- Configuration -------------
st.set_page_config(page_title="ProtHub", layout="wide")
openai.api_key = st.secrets["openai"]["api_key"]

# ------------- Constants -----------------
STRING_API_URL = "https://string-db.org/api"
STRING_OUTPUT_FORMAT = "json"
STRING_METHOD = "network"
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"

# ------------- Helper Functions -----------------
def map_sequence_to_uniprot(input_text):
    lines = input_text.strip().splitlines()
    sequence = "".join(line.strip() for line in lines if not line.startswith(">"))
    if not sequence:
        return None
    params = {"query": f'sequence:"{sequence}"', "format": "json", "size": 1}
    r = requests.get(UNIPROT_SEARCH_URL, params=params)
    if r.status_code == 200 and r.json().get("results"):
        return r.json()["results"][0]["primaryAccession"]
    return None

def get_string_interactions(uniprot_id, species=9606, min_score=0.4):
    params = {
        "identifiers": uniprot_id,
        "species": species,
        "caller_identity": "streamlit_app",
        "required_score": int(min_score * 1000),  # STRING API uses integer scores
        "output_format": STRING_OUTPUT_FORMAT,
        "method": STRING_METHOD,
    }
    r = requests.get(f"{STRING_API_URL}/{STRING_OUTPUT_FORMAT}/{STRING_METHOD}", params=params)
    if r.status_code == 200:
        return r.json()
    return None

def build_network(data):
    G = nx.Graph()
    for interaction in data:
        protein1 = interaction["preferredName_A"]
        protein2 = interaction["preferredName_B"]
        score = float(interaction["combined_score"]) / 1000  # Normalize score
        if score > 0:
            G.add_edge(protein1, protein2, weight=score)
    return G

def find_hub_genes(G, top_n=10):
    degrees = dict(G.degree())
    hub_genes = sorted(degrees, key=degrees.get, reverse=True)[:top_n]
    return hub_genes

def create_graph_figure(G, hub_genes):
    if not G.edges():
        return go.Figure(layout=go.Layout(title="No interactions found for the given parameters."))

    pos = nx.spring_layout(G, seed=42)
    degrees = dict(G.degree())

    edge_x, edge_y = [], []
    for src, dst in G.edges():
        x0, y0 = pos[src]
        x1, y1 = pos[dst]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=1.5, color='gray'),
        hoverinfo='none',
        mode='lines'
    )

    node_x, node_y, node_size, node_color, node_text, hover_text = [], [], [], [], [], []
    for node in G.nodes():
        x, y = pos[node]
        degree = degrees[node]
        node_x.append(x)
        node_y.append(y)
        node_size.append(15 + degree * 2)
        node_color.append('red' if node in hub_genes else 'royalblue')
        node_text.append(node)  # Just the node name
        hover_text.append(f"{node}<br>Degree: {degree}")

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers+text',
        text=node_text,  # Use node_text here
        textposition="middle center",
        textfont=dict(color='white', size=10),
        marker=dict(size=node_size, color=node_color, line=dict(width=1, color='white')),
        hoverinfo='text',
        hovertext=hover_text  # Use hover_text here
    )

    layout = go.Layout(
        title="Protein Interaction Network",
        titlefont=dict(color='#333'),
        showlegend=False,
        hovermode='closest',
        margin=dict(b=20, l=20, r=20, t=40),
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        plot_bgcolor='rgba(255,255,255,1)',
        paper_bgcolor='rgba(255,255,255,1)',
    )

    fig = go.Figure(data=[edge_trace, node_trace], layout=layout)
    return fig

def generate_summary(text, max_tokens=150):
    try:
        response = openai.ChatCompletion.create(
            model="gpt-3.5-turbo",
            messages=[
                {"role": "user", "content": f"Summarize: {text}"}
            ],
            max_tokens=max_tokens
        )
        return response.choices[0].message.content.strip()
    except Exception as e:
        st.error(f"Error generating summary: {e}")
        return None

# ------------- Main App -------------------
st.title("🧬 Prot'n'Hub")

input_type = st.radio("Input Type", ["UniProt ID", "Raw Sequence"])
user_input = st.text_area("Enter UniProt ID or Raw Sequence", height=100)

species_dict = {
    "Human (Homo sapiens)": 9606,
    "Mouse (Mus musculus)": 10090,
    "Rat (Rattus norvegicus)": 10116,
    "Zebrafish (Danio rerio)": 7955,
    "Fruit fly (Drosophila melanogaster)": 7227,
    "Custom (enter manually)": None,
}
selected_species = st.selectbox("Choose species", list(species_dict.keys()))
species = st.number_input("Enter NCBI Taxonomy ID:", value=9606) if selected_species == "Custom (enter manually)" else species_dict[selected_species]
score_threshold = st.slider("Minimum Interaction Score", 0.0, 1.0, 0.4, 0.05)

if st.button("Analyze Network"):
    with st.spinner("Fetching and building network..."):
        if input_type == "Raw Sequence":
            uniprot_id = map_sequence_to_uniprot(user_input.strip())
            if not uniprot_id:
                st.error("Could not map the sequence to a UniProt ID.")
                st.stop()
            else:
                st.success(f"Mapped to UniProt ID: {uniprot_id}")
        else:
            uniprot_id = user_input.strip()

        data = get_string_interactions(uniprot_id, species, score_threshold)
        if not data:
            st.error("No interaction data found.")
        else:
            G = build_network(data)
            hub_genes = find_hub_genes(G)
            st.subheader("🔗 Top Hub Genes")
            st.success(", ".join(hub_genes))

            fig = create_graph_figure(G, hub_genes)
            st.plotly_chart(fig, use_container_width=True)

            # Download button for the network (as JSON)
            network_json = nx.node_link_data(G)
            st.download_button(
                label="Download Network (JSON)",
                data=str(network_json).encode("utf-8"),
                file_name="protein_network.json",
                mime="application/json",
            )

            # Download button for the network (as PNG)
            try:
                fig.write_image("protein_network.png")
                with open("protein_network.png", "rb") as file:
                    st.download_button(
                        label="Download Network (PNG)",
                        data=file.read(),
                        file_name="protein_network.png",
                        mime="image/png",
                    )
            except Exception as e:
                st.error(f"Error saving or downloading PNG: {e}")

            # Summarize the network (using GPT)
            summary = generate_summary(str(network_json))
            if summary:
                st.subheader("📝 Network Summary")
                st.info(summary)
