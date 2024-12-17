from collections import defaultdict
import rdflib
from biocypher._logger import logger
from rdflib import Graph, RDF, URIRef, Literal
import os
from biocypher._create import BioCypherNode, BioCypherEdge

data_folder = "./"
graph = Graph()
"""Load Turtle files from the 'data' folder into the RDF graph."""
for file_name in os.listdir(data_folder):
    if file_name.endswith(".ttl"):
        file_path = os.path.join(data_folder, file_name)
        logger.info(f"Loading {file_name}")
        graph.parse(file_path, format="turtle")

def get_protein_nodes():
    protein_properties = {}

    # Define mapping from RDF URIs to human-readable property names
    property_mapping = {
        URIRef("http://www.w3.org/1999/02/22-rdf-syntax-ns#type"): "type",
        URIRef("http://purl.org/dc/terms/conformsTo"): "conforms_to",
        URIRef("http://schema.org/bioChemInteraction"): "bioChemInteraction",
        URIRef("http://schema.org/identifier"): "identifier",
        URIRef("http://schema.org/name"): "name",
        URIRef("http://schema.org/sameAs"): "sameAs",
        URIRef("http://schema.org/taxonomicRange"): "taxonomicRange",
        URIRef("http://schema.org/includedInDataset"): "included_in_dataset",
        URIRef("http://schema.org/inChIKey"): "inChIKey"
    }
    # SPARQL query to get all properties for each Protein node
    query = """
        PREFIX schema: <http://schema.org/>
        SELECT ?protein ?prop ?value
        WHERE 
        {
          ?protein a schema:Protein ;
                   ?prop ?value .
        }
        """

    # Execute the query
    results = graph.query(query)

    # Collect properties for each protein
    for row in results:
        protein = row[0]  # The protein URI
        prop = row[1]  # The property URI
        value = row[2]  # The value of the property

        # Initialize dictionary for the protein if it doesn't exist
        if protein not in protein_properties:
            protein_properties[protein] = {}

        # Map the RDF property URI to a human-readable name if it exists in the mapping
        readable_prop = property_mapping.get(prop)

        if readable_prop:
            # Initialize list for the property if not already initialized
            if readable_prop not in protein_properties[protein]:
                protein_properties[protein][readable_prop] = []

            # Convert value to string if it is Literal, leave URIRefs as-is
            if isinstance(value, Literal):
                protein_properties[protein][readable_prop].append(str(value))
            elif isinstance(value, URIRef):
                protein_properties[protein][readable_prop].append(value)


    proteins=[]
    # Create BioCypherNode objects for each protein based on its properties
    for protein, properties in protein_properties.items():
        # Set node_id and preferred_id separately, ensuring they are not in bio_properties
        protein_id = str(protein)  # `node_id`
        preferred_id = properties.get("name", [None])[0]  # First identifier, if it exists

        # Prepare additional properties, excluding 'id' and 'preferred_id'
        bio_properties = {key: values for key, values in properties.items() if key != "name"}

        # Create the BioCypherNode without id and preferred_id in properties
        bio_node = BioCypherNode(
            node_id=protein_id,
            node_label="protein",
            preferred_id=preferred_id,
            properties=bio_properties,
        )

        # Print node info for debugging, excluding `id` and `preferred_id` from properties
        print(f"Node ID: {bio_node.node_id}")
        print(f"Node Label: {bio_node.node_label}")
        print("Properties:")
        for prop, values in bio_node.properties.items():
            print(f"  Property: {prop}")
            print(f"    - {values}")
        print("\n")  # Separation for readability

        proteins.append(bio_node)
    return proteins


def get_biochemicalentity_nodes():
    """
    Returns a list of BioCypherNode objects for BioChemEntity nodes with human-readable property names.
    """
    print("Generating nodes for BioChemEntity with human-readable property names.")

    # Define mapping from RDF URIs to human-readable property names
    property_mapping = {
        URIRef("http://www.w3.org/1999/02/22-rdf-syntax-ns#type"): "type",
        URIRef("http://schema.org/name"): "name",
        URIRef("http://schema.org/identifier"): "brenda_identifier",
        URIRef("http://schema.org/includedInDataset"): "included_in_dataset",
    }

    # Dictionary to store properties of each BioChemEntity node
    biochementity_properties = {}

    # Updated SPARQL query to get properties for BioChemEntity and exclude unwanted node
    query = """
    PREFIX schema: <http://schema.org/>
    SELECT ?biochementity ?prop ?value
    WHERE 
    {
      ?biochementity a schema:BioChemEntity ;
                     ?prop ?value .
      FILTER(?biochementity != <https://bioschemas.org/profiles/MolecularEntity/0.5-RELEASE>)
    }
    """

    # Execute the query
    results = graph.query(query)

    # Collect properties for each BioChemEntity
    for row in results:
        biochementity = row[0]  # The BioChemEntity URI
        prop = row[1]  # The property URI
        value = row[2]  # The value of the property

        # Initialize dictionary for the BioChemEntity if it doesn't exist
        if biochementity not in biochementity_properties:
            biochementity_properties[biochementity] = {}

        # Map the RDF property URI to a human-readable name using property_mapping
        readable_prop = property_mapping.get(prop, str(prop))  # Default to URI if no mapping

        # Initialize list for the property if not already initialized
        if readable_prop not in biochementity_properties[biochementity]:
            biochementity_properties[biochementity][readable_prop] = []

        # Convert value to string if it's a Literal, otherwise leave URIRefs as they are
        if isinstance(value, Literal):
            biochementity_properties[biochementity][readable_prop].append(str(value))
        elif isinstance(value, URIRef):
            biochementity_properties[biochementity][readable_prop].append(value)

    # Create BioCypherNode objects for each BioChemEntity based on its properties
    biochems = []
    for biochementity, properties in biochementity_properties.items():
        # Set node_id and preferred_id separately
        biochementity_id = str(biochementity)  # `node_id`
        preferred_id = properties.get("name", [None])[0]  # Use first name as preferred_id, if it exists

        # Prepare additional properties (mapped to human-readable names)
        bio_properties = {key: values for key, values in properties.items() if key != "name"}

        # Create the BioCypherNode with human-readable property names
        bio_node = BioCypherNode(
            node_id=biochementity_id,
            node_label="biochemicalentity",
            preferred_id=preferred_id,
            properties=bio_properties,
        )

        # Print node info for debugging
        print(f"Node ID: {bio_node.node_id}")
        print(f"Node Label: {bio_node.node_label}")
        print("Properties:")
        for prop, values in bio_node.properties.items():
            print(f"  Property: {prop}")
            for value in values:
                print(f"    - {value}")
        print("\n")  # Separation for readability

        biochems.append(bio_node)

    return biochems

def get_molecularentity_nodes():
    """
    Returns a list of BioCypherNode objects for MolecularEntity nodes with human-readable property names.
    """
    print("Generating nodes for MolecularEntity with human-readable property names.")

    # Define mapping from RDF URIs to human-readable property names
    property_mapping = {
        URIRef("http://www.w3.org/1999/02/22-rdf-syntax-ns#type"): "type",
        URIRef("http://schema.org/name"): "name",
        URIRef("http://schema.org/identifier"): "identifier",
        URIRef("http://schema.org/includedInDataset"): "included_in_dataset",
        URIRef("http://schema.org/url"): "url",
        URIRef("http://schema.org/inChI"): "InChI",
        URIRef("http://schema.org/inChIKey"): "InChIKey",
        URIRef("http://schema.org/iupacName"): "iupac_name",
        URIRef("http://schema.org/smiles"): "smiles",
        URIRef("http://schema.org/molecularFormula"): "molecular_formula",
        URIRef("http://schema.org/monoisotopicMolecularWeight"): "monoisotopic_molecular_weight",
        URIRef("http://schema.org/molecularWeight"): "molecular_weight",
        URIRef("http://schema.org/image"): "image",
        URIRef("http://schema.org/description"): "description",
        URIRef("http://schema.org/hasRepresentation"): "has_representation",
        URIRef("http://purl.org/dc/terms/conformsTo"): "conforms_to",
    }

    # Dictionary to store properties of each MolecularEntity node
    molecularentity_properties = {}

    # SPARQL query to get properties for MolecularEntity, excluding unwanted profiles
    query = """
    PREFIX schema: <http://schema.org/>
    SELECT ?molecularentity ?prop ?value
    WHERE 
    {
      ?molecularentity a schema:MolecularEntity ;
                       ?prop ?value .
    }
    """

    # Execute the query
    results = graph.query(query)

    # Collect properties for each MolecularEntity
    for row in results:
        molecularentity = row[0]  # The MolecularEntity URI
        prop = row[1]  # The property URI
        value = row[2]  # The value of the property

        # Initialize dictionary for the MolecularEntity if it doesn't exist
        if molecularentity not in molecularentity_properties:
            molecularentity_properties[molecularentity] = {}

        # Map the RDF property URI to a human-readable name using property_mapping
        readable_prop = property_mapping.get(prop, str(prop))  # Default to URI if no mapping

        # Initialize list for the property if not already initialized
        if readable_prop not in molecularentity_properties[molecularentity]:
            molecularentity_properties[molecularentity][readable_prop] = []

        # Convert value to string if it's a Literal, otherwise leave URIRefs as they are
        if isinstance(value, Literal):
            molecularentity_properties[molecularentity][readable_prop].append(str(value))
        elif isinstance(value, URIRef):
            molecularentity_properties[molecularentity][readable_prop].append(str(value))

    # Create BioCypherNode objects for each MolecularEntity based on its properties
    molecularentities = []
    for molecularentity, properties in molecularentity_properties.items():
        # Set node_id and preferred_id separately (no specific preferred_id logic, use first value if available)
        molecularentity_id = str(molecularentity)  # `node_id`
        #preferred_id = properties.get("name", [None])[0] if properties.get("name") else None  # Safe access to 'name'

        # Prepare additional properties (mapped to human-readable names)
        bio_properties = {}
        for key, values in properties.items():
            # Ensure we don't get an error when iterating over None or empty lists
            if values:
                bio_properties[key] = values

        # Create the BioCypherNode with human-readable property names
        bio_node = BioCypherNode(
            node_id=molecularentity_id,
            node_label="molecularentity",
            preferred_id="name",
            properties=bio_properties,
        )

        # Print node info for debugging
        print(f"Node ID: {bio_node.node_id}")
        print(f"Node Label: {bio_node.node_label}")
        print("Properties:")
        for prop, values in bio_node.properties.items():
            print(f"  Property: {prop}")
            for value in values:
                print(f"    - {value}")
        print("\n")  # Separation for readability

        molecularentities.append(bio_node)

    return molecularentities

def get_dataset_nodes():
    dataset_properties = {}

    # Define the exclusive property mappings from the provided list
    property_mapping = {
        URIRef("http://www.w3.org/1999/02/22-rdf-syntax-ns#type"): "type",
        URIRef("http://purl.org/dc/terms/conformsTo"): "conforms_to",
        URIRef("http://schema.org/name"): "name",
        URIRef("http://schema.org/description"): "description",
        URIRef("http://schema.org/datePublished"): "datePublished",
        URIRef("http://schema.org/keywords"): "keywords",
        URIRef("http://schema.org/license"): "license",
        URIRef("http://schema.org/inLanguage"): "inLanguage",
        URIRef("http://schema.org/url"): "url",
        URIRef("http://schema.org/identifier"): "identifier",
        URIRef("http://schema.org/citation"): "citation",
        URIRef("http://schema.org/measurementTechnique"): "measurementTechnique",
        URIRef("http://schema.org/isAccessibleForFree"): "isAccessibleForFree",
        URIRef("http://schema.org/dateCreated"): "dateCreated",
        URIRef("http://schema.org/dateModified"): "dateModified",
        URIRef("http://schema.org/isPartOf"): "isPartOf",
        URIRef("http://schema.org/taxonomicRange"): "taxonomicRange"
    }

    # SPARQL query to get all properties for each Dataset node
    query = """
        PREFIX schema: <http://schema.org/>
        SELECT ?dataset ?prop ?value
        WHERE 
        {
          ?dataset a schema:Dataset ;
                   ?prop ?value .
        }
        """

    # Execute the query
    results = graph.query(query)

    # Collect properties for each dataset
    for row in results:
        dataset = row[0]  # The dataset URI
        prop = row[1]  # The property URI
        value = row[2]  # The value of the property

        # Initialize dictionary for the dataset if it doesn't exist
        if dataset not in dataset_properties:
            dataset_properties[dataset] = {}

        # Map the RDF property URI to a human-readable name if it exists in the mapping
        readable_prop = property_mapping.get(prop)

        if readable_prop:
            # Initialize list for the property if not already initialized
            if readable_prop not in dataset_properties[dataset]:
                dataset_properties[dataset][readable_prop] = []

            # Convert value to string if it is Literal, leave URIRefs as-is
            if isinstance(value, Literal):
                dataset_properties[dataset][readable_prop].append(str(value))
            elif isinstance(value, URIRef):
                dataset_properties[dataset][readable_prop].append(str(value))

    datasets = []
    # Create BioCypherNode objects for each dataset based on its properties
    for dataset, properties in dataset_properties.items():
        # Set node_id and preferred_id separately, ensuring they are not in dataset_properties
        dataset_id = str(dataset)  # `node_id`
        preferred_id = properties.get("name", [None])[0]  # First identifier, if it exists

        # Prepare additional properties, excluding 'id' and 'preferred_id'
        bio_properties = {key: values for key, values in properties.items() if key != "name"}

        # Create the BioCypherNode without id and preferred_id in properties
        bio_node = BioCypherNode(
            node_id=dataset_id,
            node_label="dataset",
            preferred_id=preferred_id,
            properties=bio_properties,
        )

        # Print node info for debugging, excluding `id` and `preferred_id` from properties
        print(f"Node ID: {bio_node.node_id}")
        print(f"Node Label: {bio_node.node_label}")
        print("Properties:")
        for prop, values in bio_node.properties.items():
            print(f"  Property: {prop}")
            print(f"    - {values}")
        print("\n")  # Separation for readability

        datasets.append(bio_node)

    return datasets

def get_organization_nodes():
    organization_properties = {}

    # Define the exclusive property mappings for Organization
    property_mapping = {
        URIRef("http://www.w3.org/1999/02/22-rdf-syntax-ns#type"): "type",
        URIRef("http://schema.org/name"): "name",
        URIRef("http://schema.org/url"): "url"
    }

    # SPARQL query to get all properties for each Organization node
    query = """
        PREFIX schema: <http://schema.org/>
        SELECT ?organization ?prop ?value
        WHERE 
        {
          ?organization a schema:Organization ;
                       ?prop ?value .
        }
        """

    # Execute the query
    results = graph.query(query)

    # Collect properties for each organization
    for row in results:
        organization = row[0]  # The organization URI
        prop = row[1]  # The property URI
        value = row[2]  # The value of the property

        # Initialize dictionary for the organization if it doesn't exist
        if organization not in organization_properties:
            organization_properties[organization] = {}

        # Map the RDF property URI to a human-readable name if it exists in the mapping
        readable_prop = property_mapping.get(prop)

        if readable_prop:
            # Initialize list for the property if not already initialized
            if readable_prop not in organization_properties[organization]:
                organization_properties[organization][readable_prop] = []

            # Convert value to string if it is Literal, leave URIRefs as-is
            if isinstance(value, Literal):
                organization_properties[organization][readable_prop].append(str(value))
            elif isinstance(value, URIRef):
                organization_properties[organization][readable_prop].append(str(value))

    organizations = []
    # Create BioCypherNode objects for each organization based on its properties
    for organization, properties in organization_properties.items():
        # Set node_id and preferred_id separately, ensuring they are not in organization_properties
        organization_id = str(organization)  # `node_id`
        preferred_id = properties.get("name", [None])[0]  # First identifier, if it exists

        # Prepare additional properties, excluding 'id' and 'preferred_id'
        bio_properties = {key: values for key, values in properties.items() if key != "name"}

        # Create the BioCypherNode without id and preferred_id in properties
        bio_node = BioCypherNode(
            node_id=organization_id,
            node_label="organization",
            preferred_id=preferred_id,
            properties=bio_properties,
        )

        # Print node info for debugging, excluding `id` and `preferred_id` from properties
        print(f"Node ID: {bio_node.node_id}")
        print(f"Node Label: {bio_node.node_label}")
        print("Properties:")
        for prop, values in bio_node.properties.items():
            print(f"  Property: {prop}")
            print(f"    - {values}")
        print("\n")  # Separation for readability

        organizations.append(bio_node)

    return organizations

def get_taxon_nodes():
    taxon_properties = {}

    # Define the exclusive property mappings for Taxon
    property_mapping = {
        URIRef("http://www.w3.org/1999/02/22-rdf-syntax-ns#type"): "type",
        URIRef("http://purl.org/dc/terms/conformsTo"): "conforms_to",
        URIRef("http://schema.org/includedInDataset"): "included_in_dataset",
        URIRef("http://schema.org/name"): "name",
        URIRef("dwc:vernacularName"): "vernacular_name",
        URIRef("http://schema.org/parentTaxon"): "parent_taxon",
        URIRef("http://schema.org/taxonRank"): "taxon_rank",
        URIRef("http://schema.org/identifier"): "identifier"
    }

    # SPARQL query to get all properties for each Taxon node
    query = """
        PREFIX schema: <http://schema.org/>
        PREFIX dwc: <http://rs.tdwg.org/dwc/terms/>
        SELECT ?taxon ?prop ?value
        WHERE 
        {
          ?taxon a schema:Taxon ;
                 ?prop ?value .
        }
        """

    # Execute the query
    results = graph.query(query)

    # Collect properties for each taxon
    for row in results:
        taxon = row[0]  # The taxon URI
        prop = row[1]  # The property URI
        value = row[2]  # The value of the property

        # Initialize dictionary for the taxon if it doesn't exist
        if taxon not in taxon_properties:
            taxon_properties[taxon] = {}

        # Map the RDF property URI to a human-readable name if it exists in the mapping
        readable_prop = property_mapping.get(prop)

        if readable_prop:
            # Initialize list for the property if not already initialized
            if readable_prop not in taxon_properties[taxon]:
                taxon_properties[taxon][readable_prop] = []

            # Convert value to string if it is Literal, leave URIRefs as-is
            if isinstance(value, Literal):
                taxon_properties[taxon][readable_prop].append(str(value))
            elif isinstance(value, URIRef):
                taxon_properties[taxon][readable_prop].append(str(value))

    taxa = []
    # Create BioCypherNode objects for each taxon based on its properties
    for taxon, properties in taxon_properties.items():
        # Set node_id and preferred_id separately, ensuring they are not in taxon_properties
        taxon_id = str(taxon)  # `node_id`
        preferred_id = properties.get("name", [None])[0]  # First identifier, if it exists

        # Prepare additional properties, excluding 'id' and 'preferred_id'
        bio_properties = {key: values for key, values in properties.items() if key != "name"}

        # Create the BioCypherNode without id and preferred_id in properties
        bio_node = BioCypherNode(
            node_id=taxon_id,
            node_label="taxon",
            preferred_id=preferred_id,
            properties=bio_properties,
        )

        # Print node info for debugging, excluding `id` and `preferred_id` from properties
        print(f"Node ID: {bio_node.node_id}")
        print(f"Node Label: {bio_node.node_label}")
        print("Properties:")
        for prop, values in bio_node.properties.items():
            print(f"  Property: {prop}")
            print(f"    - {values}")
        print("\n")  # Separation for readability

        taxa.append(bio_node)

    return taxa

def extract_properties_from_graph(graph):
    """
    Extracts all property URIs and one example value from the given RDF graph.
    This function looks for all properties of type `schema:MolecularEntity`.

    Args:
        graph (rdflib.Graph): The RDF graph to extract properties from.

    Returns:
        dict: A dictionary with property URIs as keys and list of example values (from the file).
    """
    properties = {}

    # Query to get all properties and their values for MolecularEntity nodes
    query = """
            PREFIX schema: <http://schema.org/>
            SELECT ?taxon ?prop ?value
            WHERE 
            {
                ?taxon a schema:Taxon  ;
                    ?prop ?value .
            }
            """

    results = graph.query(query)

    # Collect properties for each MolecularEntity
    for row in results:
        molecularentity = row[0]  # The MolecularEntity URI
        prop = row[1]  # The property URI
        value = row[2]  # The value of the property

        # Initialize the property list if it's not already in the dictionary
        if prop not in properties:
            properties[prop] = []

        # Store only the first example value for each property
        if not properties[prop]:  # Add only if no examples have been recorded yet
            # Store the value (either as string for Literals or URIRef for URIRefs)
            if isinstance(value, Literal):
                properties[prop].append(str(value))
            elif isinstance(value, URIRef):
                properties[prop].append(str(value))
            #else:
                #print("EEERRRROR ",prop,value)

    return properties


def extract_properties_from_ttl(ttl_file_path):
    """
    Extract properties and examples from a TTL file (Turtle format).

    Args:
        ttl_file_path (str): The file path to the TTL file.

    Returns:
        dict: A dictionary containing properties and example values (from this file).
    """
    graph = rdflib.Graph()
    graph.parse(ttl_file_path, format='ttl')  # Parse the TTL file into a graph

    # Extract properties from the RDF graph
    return extract_properties_from_graph(graph)


def summarize_properties(ttl_files):
    """
    Processes multiple TTL files, extracts properties, and displays the results.

    Args:
        ttl_files (list): List of paths to TTL files to process.

    Prints:
        The properties and their example values (one per file).
    """
    all_properties = {}

    # Process each TTL file
    for ttl_file in ttl_files:
        print(f"Processing TTL file: {ttl_file}")
        properties = extract_properties_from_ttl(ttl_file)

        # Merge properties from the current file into the global dictionary
        for prop, values in properties.items():
            if prop not in all_properties:
                all_properties[prop] = []
            # Add the value for this property from the current TTL file
            #print(values)
            if len(values) > 0:
                all_properties[prop].append(f'"{ttl_file}": {values[0]}')  # Store file-specific example

    # Print the summarized properties and examples
    print("\nSummarized properties and examples across all databases:")
    for prop, examples in all_properties.items():
        # Combine all the examples for this property across the files
        print(f"Property: {prop}")
        print(f"  Examples: {', '.join(examples)}")
        print("\n")

# Example usage:
ttl_files = [
    "coconut.ttl",  # Replace with the actual file paths
    "edal.ttl",
    "massbank.ttl",
    "metanetx.ttl",
    "nmrxiv.ttl",
    "sabio-rk.ttl",
    "wikipathways.ttl"
]

#summarize_properties(ttl_files)

get_taxon_nodes()
#get_organization_nodes()
#get_dataset_nodes()
#get_molecularentity_nodes()
#get_biochemicalentity_nodes()
#get_protein_nodes()