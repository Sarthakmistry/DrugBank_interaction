from dotenv import load_dotenv
from neo4j import GraphDatabase
from lxml import etree
import re, os

# Load .env 
load_dotenv()

URI      = os.getenv("NEO4J_URI")
USER     = os.getenv("NEO4J_USERNAME")
PASSWORD = os.getenv("NEO4J_PASSWORD")
DATABASE = os.getenv("NEO4J_DATABASE", "neo4j")

# Config
XML_PATH   = "drugbank_data.xml"
BATCH_SIZE = 200

# XML namespace setup
NS  = "http://www.drugbank.ca"
TAG = lambda t: f"{{{NS}}}{t}"

tag_drug        = TAG("drug")
tag_drugbank    = TAG("drugbank")
tag_dbid        = TAG("drugbank-id")
tag_name        = TAG("name")
tag_desc        = TAG("description")
tag_group       = TAG("group")
tag_category    = TAG("category")
tag_indication  = TAG("indication")
tag_moa         = TAG("mechanism-of-action")
tag_halflife    = TAG("half-life")
tag_cas         = TAG("cas-number")
tag_state       = TAG("state")
tag_prop        = TAG("property")
tag_kind        = TAG("kind")
tag_value       = TAG("value")
tag_iact        = TAG("drug-interaction")
tag_target_el   = TAG("target")
tag_targets_el  = TAG("targets")
tag_enzyme_el   = TAG("enzyme")
tag_enzymes_el  = TAG("enzymes")
tag_pathway_el  = TAG("pathway")
tag_pathways_el = TAG("pathways")
tag_smpdb       = TAG("smpdb-id")
tag_action      = TAG("action")
tag_polypep     = TAG("polypeptide")
tag_gene        = TAG("gene-name")
tag_ext_id      = TAG("external-identifier")
tag_resource    = TAG("resource")

# Severity classifier 
SEVERITY_RULES = [
    (re.compile(r"\b(major|severe|serious|life.?threatening|fatal|death|"
                r"hemorrhag|arrhythm|serotonin syndrome|torsade)\b", re.I), "Major"),
    (re.compile(r"\b(moderate|significant|increased risk|decreased efficacy|caution)\b", re.I), "Moderate"),
]

def classify_severity(desc):
    for pattern, label in SEVERITY_RULES:
        if pattern.search(desc or ""):
            return label
    return "Minor"

#  XML helpers 
def text(el, tag):
    found = el.find(tag)
    return (found.text or "").strip() if found is not None else ""

def calc_props(drug_el):
    formula = mol_weight = ""
    for p in drug_el.findall(f".//{tag_prop}"):
        k = text(p, tag_kind)
        v = text(p, tag_value)
        if k == "Molecular Formula": formula     = v
        if k == "Molecular Weight":  mol_weight  = v
    return formula, mol_weight

def parse_targets(drug_el, section_tag, item_tag):
    results = []
    for t in drug_el.findall(f".//{section_tag}/{item_tag}"):
        t_id    = text(t, TAG("id"))
        t_name  = text(t, tag_name)
        actions = [a.text.strip() for a in t.findall(f".//{tag_action}") if a.text]
        gene, uniprot = "", ""
        poly = t.find(f".//{tag_polypep}")
        if poly is not None:
            gene = text(poly, tag_gene)
            for ext in poly.findall(f".//{tag_ext_id}"):
                if text(ext, tag_resource) == "UniProtKB":
                    uniprot = text(ext, TAG("identifier"))
                    break
        if t_id:
            results.append({
                "id": t_id, "name": t_name,
                "gene": gene, "uniprot": uniprot,
                "actions": "|".join(actions),
            })
    return results

# Cypher queries
MERGE_DRUG = """
UNWIND $batch AS d
MERGE (n:Drug {id: d.id})
SET n.name        = d.name,
    n.description = d.description,
    n.drugType    = d.drugType,
    n.categories  = d.categories,
    n.formula     = d.formula,
    n.molWeight   = toFloat(d.molWeight),
    n.indication  = d.indication,
    n.mechanism   = d.mechanism,
    n.halfLife    = d.halfLife,
    n.casNumber   = d.casNumber,
    n.state       = d.state,
    n.approved    = true
"""

MERGE_INTERACTIONS = """
UNWIND $batch AS r
MATCH (a:Drug {id: r.src})
MATCH (b:Drug {id: r.tgt})
MERGE (a)-[rel:INTERACTS_WITH]->(b)
SET rel.severity    = r.severity,
    rel.description = r.description
"""

MERGE_TARGETS = """
UNWIND $batch AS t
MERGE (p:ProteinTarget {id: t.id})
SET p.name      = t.name,
    p.gene      = t.gene,
    p.uniprotId = t.uniprot
WITH p, t
MATCH (d:Drug {id: t.drug_id})
MERGE (d)-[r:TARGETS]->(p)
SET r.actions = t.actions
"""

MERGE_ENZYMES = """
UNWIND $batch AS t
MERGE (p:ProteinTarget {id: t.id})
SET p.name      = t.name,
    p.gene      = t.gene,
    p.uniprotId = t.uniprot
WITH p, t
MATCH (d:Drug {id: t.drug_id})
MERGE (d)-[r:METABOLIZED_BY]->(p)
SET r.actions = t.actions
"""

MERGE_PATHWAYS = """
UNWIND $batch AS pw
MERGE (p:Pathway {id: pw.id})
SET p.name     = pw.name,
    p.category = pw.category
WITH p, pw
MATCH (d:Drug {id: pw.drug_id})
MERGE (d)-[:PART_OF]->(p)
"""

CONSTRAINTS = [
    "CREATE CONSTRAINT drug_id    IF NOT EXISTS FOR (d:Drug)          REQUIRE d.id IS UNIQUE",
    "CREATE CONSTRAINT target_id  IF NOT EXISTS FOR (t:ProteinTarget) REQUIRE t.id IS UNIQUE",
    "CREATE CONSTRAINT pathway_id IF NOT EXISTS FOR (p:Pathway)       REQUIRE p.id IS UNIQUE",
    "CREATE INDEX drug_name       IF NOT EXISTS FOR (d:Drug)          ON (d.name)",
    "CREATE INDEX drug_approved   IF NOT EXISTS FOR (d:Drug)          ON (d.approved)",
    "CREATE INDEX iact_severity   IF NOT EXISTS FOR ()-[r:INTERACTS_WITH]-() ON (r.severity)",
]

# ── Loader ────────────────────────────────────────────────────────────────────
class Loader:
    def __init__(self):
        self.driver = GraphDatabase.driver(URI, auth=(USER, PASSWORD))
        print(f"Driver created → {URI}")

    def verify(self):
        self.driver.verify_connectivity()
        print("Connection verified")

    def close(self):
        self.driver.close()
        print("Driver closed.")

    def run_setup(self):
        with self.driver.session(database=DATABASE) as s:
            for q in CONSTRAINTS:
                s.run(q)
        print(" Constraints and indexes ready.")

    def flush(self, cypher, batch, label):
        if not batch:
            return
        with self.driver.session(database=DATABASE) as s:
            s.run(cypher, batch=batch)
        print(f"-> wrote {len(batch):>4} {label}")

# ── Main ingest ───────────────────────────────────────────────────────────────
def run():
    loader = Loader()
    loader.verify()
    loader.run_setup()

    drug_batch    = []
    iact_batch    = []
    target_batch  = []
    enzyme_batch  = []
    pathway_batch = []
    seen_edges    = set()
    total = {"drugs": 0, "interactions": 0, "targets": 0, "enzymes": 0, "pathways": 0}

    ctx = etree.iterparse(XML_PATH, events=("end",), tag=tag_drug)

    for _, el in ctx:
        parent = el.getparent()
        if parent is None or parent.tag != tag_drugbank:
            el.clear()
            continue

        # Filter: approved only
        groups = [g.text.strip() for g in el.findall(f".//{tag_group}") if g.text]
        if "approved" not in groups:
            el.clear()
            continue

        # Primary DrugBank ID
        drug_id = None
        for id_el in el.findall(tag_dbid):
            if id_el.get("primary") == "true":
                drug_id = (id_el.text or "").strip()
                break
        if not drug_id:
            el.clear()
            continue

        # Drug node
        formula, mol_weight = calc_props(el)
        categories = [
            c.text.strip()
            for c in el.findall(f".//{tag_category}/{tag_category}")
            if c.text
        ]
        drug_batch.append({
            "id":          drug_id,
            "name":        text(el, tag_name),
            "description": text(el, tag_desc)[:500],
            "drugType":    el.get("type", ""),
            "categories":  "|".join(categories[:3]),
            "formula":     formula,
            "molWeight":   mol_weight,
            "indication":  text(el, tag_indication)[:300],
            "mechanism":   text(el, tag_moa)[:300],
            "halfLife":    text(el, tag_halflife),
            "casNumber":   text(el, tag_cas),
            "state":       text(el, tag_state),
        })
        total["drugs"] += 1

        # Interactions — Major/Moderate only
        for iact in el.findall(f".//{tag_iact}"):
            tgt_id   = text(iact, tag_dbid)
            if not tgt_id:
                continue
            desc     = text(iact, tag_desc)
            severity = classify_severity(desc)
            if severity not in ("Major", "Moderate"):
                continue
            edge_key = tuple(sorted([drug_id, tgt_id]))
            if edge_key in seen_edges:
                continue
            seen_edges.add(edge_key)
            iact_batch.append({
                "src": drug_id, "tgt": tgt_id,
                "severity": severity, "description": desc[:500],
            })
            total["interactions"] += 1

        # Targets
        for t in parse_targets(el, tag_targets_el, tag_target_el):
            t["drug_id"] = drug_id
            target_batch.append(t)
            total["targets"] += 1

        # Enzymes
        for e in parse_targets(el, tag_enzymes_el, tag_enzyme_el):
            e["drug_id"] = drug_id
            enzyme_batch.append(e)
            total["enzymes"] += 1

        # Pathways
        for pw in el.findall(f".//{tag_pathways_el}/{tag_pathway_el}"):
            pw_id = text(pw, tag_smpdb)
            if pw_id:
                pathway_batch.append({
                    "id":       pw_id,
                    "name":     text(pw, tag_name),
                    "category": text(pw, TAG("category")),
                    "drug_id":  drug_id,
                })
                total["pathways"] += 1

        # Flush full batches
        if len(drug_batch)    >= BATCH_SIZE: loader.flush(MERGE_DRUG,          drug_batch,    "drugs");    drug_batch.clear()
        if len(iact_batch)    >= BATCH_SIZE: loader.flush(MERGE_INTERACTIONS,  iact_batch,    "interactions"); iact_batch.clear()
        if len(target_batch)  >= BATCH_SIZE: loader.flush(MERGE_TARGETS,       target_batch,  "targets");  target_batch.clear()
        if len(enzyme_batch)  >= BATCH_SIZE: loader.flush(MERGE_ENZYMES,       enzyme_batch,  "enzymes");  enzyme_batch.clear()
        if len(pathway_batch) >= BATCH_SIZE: loader.flush(MERGE_PATHWAYS,      pathway_batch, "pathways"); pathway_batch.clear()

        # Free memory
        el.clear()
        while el.getprevious() is not None:
            del el.getparent()[0]

        if total["drugs"] % 200 == 0:
            print(f"  Progress: {total}")

    # Final flush
    loader.flush(MERGE_DRUG,         drug_batch,    "drugs (final)")
    loader.flush(MERGE_INTERACTIONS, iact_batch,    "interactions (final)")
    loader.flush(MERGE_TARGETS,      target_batch,  "targets (final)")
    loader.flush(MERGE_ENZYMES,      enzyme_batch,  "enzymes (final)")
    loader.flush(MERGE_PATHWAYS,     pathway_batch, "pathways (final)")

    loader.close()
    print(f"\n Ingest complete: {total}")

# ── Run ───────────────────────────────────────────────────────────────────────
run()