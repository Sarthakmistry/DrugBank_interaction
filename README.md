# DDI Explorer

A full-stack web application for exploring drug-drug interactions using a Neo4j graph database with DrugBank data.

<img width="1849" height="910" alt="Image" src="https://github.com/user-attachments/assets/46c2a2d1-7591-4525-bf5b-654d8c1a3c64" />

---

## Stack

- **Backend:** Node.js + Express + Neo4j JavaScript Driver
- **Frontend:** React 18 (Vite) + TailwindCSS + D3.js + Recharts
- **Database:** Neo4j AuraDB (free tier) — pre-loaded with DrugBank data

---

## Database

The Neo4j graph contains filtered DrugBank data (approved drugs with Major/Moderate interactions only):

| Label | Count | Key Properties |
|---|---|---|
| `Drug` | ~4,795 | `id`, `name`, `description`, `drugType`, `categories`, `formula`, `indication`, `mechanism` |
| `ProteinTarget` | ~3,272 | `id`, `name`, `gene`, `uniprotId` |
| `Pathway` | ~48,622 | `id`, `name`, `category` |

| Relationship | Description |
|---|---|
| `INTERACTS_WITH` | Drug ↔ Drug, with `severity` (Major/Moderate) and `description` |
| `TARGETS` | Drug → ProteinTarget |
| `METABOLIZED_BY` | Drug → ProteinTarget |
| `PART_OF` | Drug → Pathway |

---

## Features

- **Drug Search** — debounced autocomplete with full drug detail panel (properties, interactions, targets, pathways)
- **Interaction Checker** — build a regimen of up to 10 drugs and check all pairwise interactions
- **Graph Explorer** — D3.js force-directed graph showing a drug's interaction network (1-hop or 2-hop), with zoom/pan and click-to-explore
- **CRUD Manager** — create, edit, and delete drugs and interactions.
- **Analytics Dashboard** — summary stats, top targeted proteins, most interactive drugs, and severity breakdown charts

---

## Project Structure

```
.
├── ingest.py           # DrugBank XML parser and Neo4j loader
├── explore.ipynb       # Cypher query exploration notebook
├── .env.example        # Environment variable template
└── ddi-explorer/
    ├── package.json    # Root — runs server + client concurrently
    ├── server/
    │   ├── index.js    # Express app
    │   ├── db.js       # Neo4j driver singleton
    │   └── routes/
    │       ├── drugs.js
    │       ├── interactions.js
    │       ├── graph.js
    │       ├── targets.js
    │       ├── pathways.js
    │       └── analytics.js
    └── client/
        ├── vite.config.js
        └── src/
            ├── App.jsx
            ├── api.js
            └── pages/
                ├── DrugSearch.jsx
                ├── InteractionChecker.jsx
                ├── GraphExplorer.jsx
                ├── CRUDManager.jsx
                └── Analytics.jsx
```

---

## Setup

### 1. Environment

Copy `.env.example` to `.env` and fill in your Neo4j AuraDB credentials:

```
NEO4J_URI=neo4j+s://<instance-id>.databases.neo4j.io
NEO4J_USERNAME=<username>
NEO4J_PASSWORD=<password>
NEO4J_DATABASE=<database-name>
PORT=3001
```

### 2. Install dependencies

```bash
cd ddi-explorer
npm install
cd client && npm install
```

### 3. Run

```bash
cd ddi-explorer
npm run dev
```

Opens the API on `http://localhost:3001` and the UI on `http://localhost:5173`.

---

## Notes

- AuraDB free-tier instances pause after ~3 days of inactivity. Resume at [console.neo4j.io](https://console.neo4j.io) before running.
- Data is filtered to approved drugs only. The DrugBank XML source file is not included in this repo.

---

## Citation

Wishart DS, Knox C, Guo AC, Cheng D, Shrivastava S, Tzur D, Gautam B, Hassanali M. DrugBank: a knowledgebase for drugs, drug actions and drug targets. *Nucleic Acids Res.* 2008 Jan;36(Database issue):D901-6.
