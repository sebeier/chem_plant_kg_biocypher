import re
from typing import Optional
from biocypher._logger import logger
from rdflib import Graph, URIRef, Literal
import os
from biocypher._create import BioCypherNode, BioCypherEdge
import hashlib

logger.debug(f"Loading module {__name__}.")

class TestAdapter:
    """
    Test BioCypher adapter. Generates nodes and edges for creating a
    knowledge graph.

    Args:
        node_types: List of node types to include in the result.
        node_fields: List of node fields to include in the result.
        edge_types: List of edge types to include in the result.
        edge_fields: List of edge fields to include in the result.
    """
    data_folder = "data"

    def __init__(
        self,
        node_types: Optional[list] = None,
        node_fields: Optional[list] = None,
        edge_types: Optional[list] = None,
        edge_fields: Optional[list] = None,
    ):
        self.graph = Graph()
        self.load_data()
        self.nodes=[node for node in self.prepare_nodes()]
        #self._set_types_and_fields(node_types, node_fields, edge_types, edge_fields)

    def load_data(self):
        """Load Turtle files from the 'data' folder into the RDF graph."""
        for file_name in os.listdir(self.data_folder):
            if file_name.endswith(".ttl"):
                file_path = os.path.join(self.data_folder, file_name)
                logger.info(f"Loading {file_name}")
                self.graph.parse(file_path, format="turtle")

    def get_taxon_nodes(self):
        print("Generating nodes for Taxon with human-readable property names.")
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
        results = self.graph.query(query)

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
            #print(f"Node ID: {bio_node.node_id}")
            #print(f"Node Label: {bio_node.node_label}")
            #print("Properties:")
            #for prop, values in bio_node.properties.items():
                #print(f"  Property: {prop}")
                #print(f"    - {values}")
            #print("\n")  # Separation for readability

            taxa.append(bio_node)

        return taxa

    def get_dataset_nodes(self):
        print("Generating nodes for Dataset with human-readable property names.")
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
        results = self.graph.query(query)

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
            #print(f"Node ID: {bio_node.node_id}")
            #print(f"Node Label: {bio_node.node_label}")
            #print("Properties:")
            #for prop, values in bio_node.properties.items():
                #print(f"  Property: {prop}")
                #print(f"    - {values}")
            #print("\n")  # Separation for readability

            datasets.append(bio_node)

        return datasets

    def get_molecular_entity_nodes(self):
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
        molecular_entity_properties = {}

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
        results = self.graph.query(query)

        # Collect properties for each MolecularEntity
        for row in results:
            molecular_entity = row[0]  # The MolecularEntity URI
            prop = row[1]  # The property URI
            value = row[2]  # The value of the property

            # Initialize dictionary for the MolecularEntity if it doesn't exist
            if molecular_entity not in molecular_entity_properties:
                molecular_entity_properties[molecular_entity] = {}

            # Map the RDF property URI to a human-readable name using property_mapping
            readable_prop = property_mapping.get(prop, str(prop))  # Default to URI if no mapping

            # Initialize list for the property if not already initialized
            if readable_prop not in molecular_entity_properties[molecular_entity]:
                molecular_entity_properties[molecular_entity][readable_prop] = []

            # Convert value to string if it's a Literal, otherwise leave URIRefs as they are
            if isinstance(value, Literal):
                molecular_entity_properties[molecular_entity][readable_prop].append(str(value))
            elif isinstance(value, URIRef):
                molecular_entity_properties[molecular_entity][readable_prop].append(str(value))

        # Create BioCypherNode objects for each MolecularEntity based on its properties
        molecular_entities = []
        for molecular_entity, properties in molecular_entity_properties.items():
            # Set node_id and preferred_id separately (no specific preferred_id logic, use first value if available)
            molecular_entity_id = str(molecular_entity)  # `node_id`
            # preferred_id = properties.get("name", [None])[0] if properties.get("name") else None  # Safe access to 'name'

            # Prepare additional properties (mapped to human-readable names)
            bio_properties = {}
            for key, values in properties.items():
                # Ensure we don't get an error when iterating over None or empty lists
                if values:
                    bio_properties[key] = values

            # Create the BioCypherNode with human-readable property names
            bio_node = BioCypherNode(
                node_id=molecular_entity_id,
                node_label="molecularentity",
                preferred_id="name",
                properties=bio_properties,
            )
            # Print node info for debugging
            #print(f"Node ID: {bio_node.node_id}")
           # print(f"Node Label: {bio_node.node_label}")
            #print("Properties:")
            #for prop, values in bio_node.properties.items():
             #   print(f"  Property: {prop}")
              #  for value in values:
               #     print(f"    - {value}")
            #print("\n")  # Separation for readability

            molecular_entities.append(bio_node)

        return molecular_entities

    def get_biochemical_entity_nodes(self):
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
        biochemical_entity_properties = {}

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
        results = self.graph.query(query)

        # Collect properties for each BioChemEntity
        for row in results:
            biochemical_entity = row[0]  # The BioChemEntity URI
            prop = row[1]  # The property URI
            value = row[2]  # The value of the property

            # Initialize dictionary for the BioChemEntity if it doesn't exist
            if biochemical_entity not in biochemical_entity_properties:
                biochemical_entity_properties[biochemical_entity] = {}

            # Map the RDF property URI to a human-readable name using property_mapping
            readable_prop = property_mapping.get(prop, str(prop))  # Default to URI if no mapping

            # Initialize list for the property if not already initialized
            if readable_prop not in biochemical_entity_properties[biochemical_entity]:
                biochemical_entity_properties[biochemical_entity][readable_prop] = []

            # Convert value to string if it's a Literal, otherwise leave URIRefs as they are
            if isinstance(value, Literal):
                biochemical_entity_properties[biochemical_entity][readable_prop].append(str(value))
            elif isinstance(value, URIRef):
                biochemical_entity_properties[biochemical_entity][readable_prop].append(value)

        # Create BioCypherNode objects for each BioChemEntity based on its properties
        biochemicals = []
        for biochemical_entity, properties in biochemical_entity_properties.items():
            # Set node_id and preferred_id separately
            biochemical_entity_id = str(biochemical_entity)  # `node_id`
            preferred_id = properties.get("name", [None])[0]  # Use first name as preferred_id, if it exists

            # Prepare additional properties (mapped to human-readable names)
            bio_properties = {key: values for key, values in properties.items() if key != "name"}

            # Create the BioCypherNode with human-readable property names
            bio_node = BioCypherNode(
                node_id=biochemical_entity_id,
                node_label="biochemicalentity",
                preferred_id=preferred_id,
                properties=bio_properties,
            )

            # Print node info for debugging
            #print(f"Node ID: {bio_node.node_id}")
            #print(f"Node Label: {bio_node.node_label}")
            #print("Properties:")
            #for prop, values in bio_node.properties.items():
                #print(f"  Property: {prop}")
                #for value in values:
                    #print(f"    - {value}")
            #print("\n")  # Separation for readability

            biochemicals.append(bio_node)

        return biochemicals

    def get_protein_nodes(self):
        print("Generating nodes for Protein with human-readable property names.")
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
        results = self.graph.query(query)

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

        proteins = []
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
            #print(f"Node ID: {bio_node.node_id}")
            #print(f"Node Label: {bio_node.node_label}")
            #print("Properties:")
            #for prop, values in bio_node.properties.items():
                #print(f"  Property: {prop}")
                #print(f"    - {values}")
            #print("\n")  # Separation for readability

            proteins.append(bio_node)
        return proteins

    def prepare_nodes(self):
        for protein in self.get_protein_nodes():
            yield protein
        for biochemical_entity in self.get_biochemical_entity_nodes():
            yield biochemical_entity
        for molecular_entity in self.get_molecular_entity_nodes():
            yield molecular_entity
        for dataset in self.get_dataset_nodes():
            yield dataset
        for taxon in self.get_taxon_nodes():
            yield taxon

    def get_nodes(self):
        for node in self.nodes:
            yield node

    #def populate_nodes(self):
        #self.nodes = list(self.get_nodes())

    def extract_numeric_id(self,uri):
        # Extract numeric part from a URI (e.g., prefix:12345 -> 12345)
        match = re.search(r'(\d+)$', uri)
        return match.group(1) if match else None

    def get_edges(self):
        """
        Returns a generator of edge tuples for edge types specified in the
        adapter constructor.

        Args:
            probability: Probability of generating an edge between two nodes.
        """

        logger.info("Generating edges.")
        #self.populate_nodes()
        if not self.nodes:
            raise ValueError("No nodes found. Please run get_nodes() first.")
        edges = []

        # Extract taxon identifiers in a dictionary where the key is the numeric ID and the value is the node_id
        taxon_nodes = {}
        for node in self.nodes:
            if node.node_label == "taxon" and "identifier" in node.properties:
                # Iterate through all identifiers for the taxon node
                for identifier in node.properties["identifier"]:
                    taxon_numeric_id = self.extract_numeric_id(str(identifier)) if identifier else None
                    if taxon_numeric_id:
                        # Store each identifier's numeric ID as key and node_id as value
                        taxon_nodes[taxon_numeric_id] = node.node_id

        # Iterate through protein nodes and check if their taxonomicRange matches any taxon identifier
        for protein in self.nodes:
            if protein.node_label == "protein" and "taxonomicRange" in protein.properties:
                # Iterate over all taxonomicRange values for the protein node
                for taxonomic_range_uri in protein.properties["taxonomicRange"]:
                    # Extract numeric ID from the taxonomicRange URI
                    taxon_numeric_id = self.extract_numeric_id(str(taxonomic_range_uri)) if taxonomic_range_uri else None

                    # Check if the taxon numeric ID exists in the taxon_nodes dictionary
                    if taxon_numeric_id and taxon_numeric_id in taxon_nodes:
                        # Match found, create edge
                        source_id = protein.node_id
                        target_id = taxon_nodes[taxon_numeric_id]
                        relationship_label = "has_taxonomic_range"

                        # Generate edge properties (if needed)
                        properties = {}

                        # Generate relationship ID using md5 hash
                        relationship_id = hashlib.md5(
                            f"{source_id}-{target_id}-{relationship_label}".encode('utf-8')).hexdigest()

                        # Debugging: print the edge being created
                        #print(
                            #f"Edge: {relationship_id} - Source: {source_id}, Target: {target_id}, Label: {relationship_label}")

                        # Create BioCypherEdge object
                        edge = BioCypherEdge(
                            relationship_id=relationship_id,
                            source_id=source_id,
                            target_id=target_id,
                            relationship_label=relationship_label,
                            properties=properties
                        )

                        # Yield the edge
                        yield edge

    def get_node_count(self):
        """
        Returns the number of nodes generated by the adapter.
        """
        return len(list(self.get_nodes()))
