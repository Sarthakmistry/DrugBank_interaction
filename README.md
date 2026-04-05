# Drug - Drug interaction
This application can be used for interacting with drug nodes from drugbank.


### Progress till now:

The data has been ingested in Neo4j Aura instance. Due to the limitations of free tier, Data has been filtered to include only approved drugs with major interactions.

```
- explore.ipynb # explores data and runs cypher queries post ingestion
- ingest.py     # Ingestion script used to parse the original xml data and model into nodes and edges.
- env.example   # Example environment file for savin credentials
```

## Citation

Wishart DS, Knox C, Guo AC, Cheng D, Shrivastava S, Tzur D, Gautam B, Hassanali M. DrugBank: a knowledgebase for drugs, drug actions and drug targets. Nucleic Acids Res. 2008 Jan;36(Database issue):D901-6.
