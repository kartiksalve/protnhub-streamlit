import streamlit as st
import requests
import networkx as nx
import plotly.graph_objects as go
# import base64 # No longer needed for background
# import os # Not used
import io

# -------------- Background Styling (Removed background image, kept container styling) --------------
def style_app_container():
    css = f"""
    <style>
    .block-container {{
        background-color: rgba(0, 0, 0, 0.6); /* Kept for content readability if desired */
        /* You can adjust or remove this background-color if you want a plain white background for the container */
        /* background-color: white; */ /* Example for a white container */
        padding: 2rem 3rem;
        border-radius: 1rem;
        /* backdrop-filter: blur(8px); */ /* Removed as it was tied to the background image effect */
        /* -webkit-backdrop-filter: blur(8px); */ /* Removed */
        color: #333; /* Changed to a more standard text color for better readability without dark background */
        /* color: #f0f0f0; */ /* Original color for dark background */
        max-width: 1000px;
        margin: auto;
    }}
    /* Optional: Style for the main app body if you want a specific color */
    .stApp {{
        /* background-color: #f0f2f6; */ /* Example: Light grey background for the whole app */
    }}
    </style>
    """
    st.markdown(css, unsafe_allow_html=True)

style_app_container() # Apply container styling

# -------------- Constants --------------
STRING_API_URL = "https://string-db.org/api"
STRING_OUTPUT_FORMAT = "json"
STRING_METHOD = "network"

# -------------- Functions --------------
def get_string_interactions(uniprot_id, species=9606, min_score=0.4):
    params = {
        "identifiers": uniprot_id,
        "species": species,
        "caller_identity": "streamlit_app",
        "required_score": int(min_score * 1000)
    }
    url = f"{STRING_API_URL}/{STRING_OUTPUT_FORMAT}/{STRING_METHOD}"
    response = requests.post(url, data=params)
    if response.status_code == 200:
        return response.json()
    return None

def build_network(data):
    G = nx.DiGraph()
    for interaction in data:
        p1 = interaction['preferredName_A']
        p2 = interaction['preferredName_B']
        score = interaction['score']
        G.add_edge(p1, p2, weight=score)
    return G

def find_hub_genes(G, top_n=5):
    degree_dict = dict(G.degree())
    sorted_genes = sorted(degree_dict.items(), key=lambda x: x[1], reverse=True)
    return [gene for gene, _ in sorted_genes[:top_n]]

def create_graph_figure(G, hub_genes):
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

    node_x, node_y, node_size, node_color, node_text = [], [], [], [], []
    for node in G.nodes():
        x, y = pos[node]
        degree = degrees[node]
        node_x.append(x)
        node_y.append(y)
        node_size.append(15 + degree * 2)
        node_color.append('red' if node in hub_genes else 'royalblue')
        node_text.append(f"{node}<br>Degree: {degree}")

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers+text',
        text=[node for node in G.nodes()],
        textposition="middle center",
        textfont=dict(color='white', size=10), # Text on nodes
        marker=dict(size=node_size, color=node_color, line=dict(width=1, color='white')),
        hoverinfo='text',
        hovertext=node_text
    )

    layout = go.Layout(
        title="Protein Interaction Network",
        titlefont=dict(color='#333'), # Title color
        showlegend=False,
        hovermode='closest',
        margin=dict(b=20, l=20, r=20, t=40),
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        plot_bgcolor='rgba(255,255,255,1)', # White plot background
        paper_bgcolor='rgba(255,255,255,1)', # White paper background
    )

    return go.Figure(data=[edge_trace, node_trace], layout=layout)

# -------------- Streamlit UI --------------
st.title("🧬 Prot'n'Hub – Protein Interaction & Hub Gene Explorer")

tabs = st.tabs(["Home", "About"])
with tabs[0]:
    st.header("Explore Protein Network")

    user_input = st.text_area("Enter Protein Name or UniProt ID", height=120)

    st.subheader("Species Selection")
    species_dict = {
        "Human (Homo sapiens)": 9606,
        "Mouse (Mus musculus)": 10090,
        "Rat (Rattus norvegicus)": 10116,
        "Zebrafish (Danio rerio)": 7955,
        "Fruit fly (Drosophila melanogaster)": 7227,
        "Custom (enter manually)": None,
    }
    selected_species = st.selectbox("Choose species", list(species_dict.keys()), index=0)
    if selected_species == "Custom (enter manually)":
        species = st.number_input("Enter NCBI Taxonomy ID", value=9606)
    else:
        species = species_dict[selected_species]

    score_threshold = st.slider("Interaction Score Threshold", 0.0, 1.0, 0.4, 0.05)

    if st.button("🔍 Analyze"):
        with st.spinner("Fetching and analyzing data..."):
            uniprot_id = user_input.strip()
            if not uniprot_id:
                st.warning("Please enter a Protein Name or UniProt ID.")
                st.stop()

            string_data = get_string_interactions(uniprot_id, species, score_threshold)
            if not string_data:
                st.error("❌ No interaction data found or error fetching data. Check your input and try again.")
                st.stop()

            G = build_network(string_data)
            if G.number_of_nodes() == 0:
                st.warning("ℹ️ No network could be built with the given parameters. Try a lower score threshold or a different protein.")
                st.stop()

            hub_genes = find_hub_genes(G)
            if hub_genes:
                st.success(f"Top Hub Genes: {', '.join(hub_genes)}")
            else:
                st.info("No distinct hub genes found based on the current network.")


            fig = create_graph_figure(G, hub_genes)
            st.plotly_chart(fig, use_container_width=True)

            # Save to PNG and offer download
            buf = io.BytesIO()
            fig.write_image(buf, format="png", width=1000, height=800, engine="kaleido")
            st.download_button(
                label="📥 Download Network as PNG",
                data=buf.getvalue(),
                file_name=f"{uniprot_id}_network.png",
                mime="image/png"
            )

            with st.expander("📊 Network Analysis Results", expanded=True):
                st.write(f"⭐ **Nodes**: {G.number_of_nodes()}")
                st.write(f"List of Nodes: {', '.join(list(G.nodes()))}")
                st.write(f"⭐ **Edges**: {G.number_of_edges()}")
                # st.write(f"List of Edges: {list(G.edges())}") # Can be very long
                degree_dict = dict(G.degree())
                st.write("⭐ **Node Degrees:**")
                for node, degree in sorted(degree_dict.items(), key=lambda item: item[1], reverse=True)[:10]: # Show top 10
                    st.write(f"- {node}: {degree}")
                if len(degree_dict) > 10:
                    st.write("... and more.")

                if degree_dict:
                    main_hub = max(degree_dict, key=degree_dict.get)
                    st.info(f"🏆 **Main Hub Gene (Highest Degree):** {main_hub} (Degree: {degree_dict[main_hub]})")
                else:
                    st.info("No nodes to determine a main hub gene.")


with tabs[1]:
    st.header("About Prot'n'Hub")
    st.markdown("""
    Prot'n'Hub is a Streamlit-based interactive application for exploring protein-protein interaction networks.

    **Features:**
    - Input UniProt ID or protein name
    - Visualize interaction graphs
    - Detect top hub genes
    - Interactive, styled Plotly graphs
    - Customizable species and score thresholds
    - Downloadable graph as PNG image

    **Powered By:**
    - STRING API
    - UniProt API (Implicitly via STRING for name resolution)
    - NetworkX, Plotly, and Streamlit

    ---

    ### 🧑‍🏫 Quick Guide

      Prot'n'Hub is designed to be simple and informative. Here's a quick guide:

    **🔹 Explore Protein Network:**  
    Enter a *protein name* (like "TP53") or a *UniProt ID* (like "P04637"). The app uses this to search for protein-protein interactions from the STRING database.

    **🔹 Species Selection:**  
    Choose the species your protein belongs to. For example:
    - Human = Homo sapiens
    - Mouse = Mus musculus  
    If your species isn't listed, choose **Custom** and enter its **NCBI Taxonomy ID** (a unique number for each species).

    **🔹 Interaction Score Threshold:**  
    This sets the minimum confidence score for the interactions.  
    - Lower values (e.g. 0.2) = more connections, but lower reliability  
    - Higher values (e.g. 0.7+) = fewer connections, but higher reliability  
    Default is **0.4**, which balances both.

    **🔹 Nodes:**  
    Each circle in the graph is a *protein*. The number of nodes shows how many proteins are in your interaction network.

    **🔹 Edges:**  
    Lines connecting the nodes. Each edge represents an interaction between two proteins.

    **🔹 Node Degrees:**  
    This shows how many connections (edges) each protein (node) has.  
    Proteins with high degree values are often *hub genes*—key proteins that interact with many others.

    **🔹 Main Hub Gene:**  
    The protein with the most connections in your network. These are often biologically important and can be potential targets for further research.

    **🔹 Download Graph as PNG:**  
    Click this button to save the visual interaction network as a PNG image for reports or presentations.

    ---

    ### 📄 Acknowledgement
      I would like to express my sincere gratitude to the following resources and individuals who contributed to the development of this application:

    **STRING Database**: For providing comprehensive protein-protein interaction data, which forms the core of this application's functionality.  
    **UniProt Knowledgebase**: For enabling the mapping of protein sequences to UniProt IDs, a crucial step in processing user-provided sequence input.  
    **NetworkX**: For the powerful tools used in network analysis and manipulation, allowing for the construction and processing of protein interaction graphs.  
    **Plotly**: For the creation of interactive and visually appealing network visualizations, enhancing the user experience.  
    **Streamlit**: For providing a user-friendly framework for building and deploying the web application.  

    I extend my deepest gratitude to **Dr. Kushagra Kashyap**, for his invaluable guidance, support, and expertise throughout this project. His insights and encouragement were instrumental  in shaping the application and overcoming the challenges encountered during its development.

    ---

    **Developer Information:**

    Hi, I'm Kartik Parag Salve, the creator of Prot'n'Hub.
    As someone currently pursuing my Master's in Bioinformatics at Deccan Education Society Pune University,
    I've always been fascinated by the power of protein-protein interaction networks.
    That's why I built Prot'n'Hub which is a Streamlit app that makes exploring these interactions intuitive and engaging.

    I'm a big fan of the rich data available through the STRING and UniProt APIs,
    and the amazing visualization capabilities of NetworkX and Plotly.
    Prot'n'Hub lets you:
    - Input UniProt IDs or protein names
    - Visualize interaction graphs
    - Find key hub genes
    - Customize species and score thresholds
    - Download the graph for reports or presentations

    It's been exciting to bring these tools together in a user-friendly application,
    and I hope Prot'n'Hub proves useful for others exploring the fascinating world of protein interactions.

    **Contact:** 3522411023@despu.edu.in
    """)