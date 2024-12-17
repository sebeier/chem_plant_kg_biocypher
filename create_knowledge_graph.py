from biocypher import BioCypher, Resource
from template_package.adapters.test_adapter import (
    TestAdapter
   #ExampleAdapterNodeType,
    #ExampleAdapterEdgeType,
    #ExampleAdapterProteinField,
    #ExampleAdapterDiseaseField,
)

# Instantiate the BioCypher interface
# You can use `config/biocypher_config.yaml` to configure the framework or
# supply settings via parameters below
bc = BioCypher()

# Download and cache resources (change the directory in the options if needed)


# Create a protein adapter instance
adapter = TestAdapter(

    # we can leave edge fields empty, defaulting to all fields in the adapter
)

# Create a knowledge graph from the adapter
bc.write_nodes(adapter.get_nodes())
bc.write_edges(adapter.get_edges())

# Write admin import statement
bc.write_import_call()

# Print summary
bc.summary()
