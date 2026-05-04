import { Router } from 'express';
import driver, { session as openSession } from '../db.js';

const router = Router();

// GET /api/drugs/search?q=<name>
router.get('/search', async (req, res) => {
  const { q } = req.query;
  if (!q || q.length < 2) return res.json([]);

  const session = openSession();
  try {
    const result = await session.run(
      `MATCH (d:Drug)
       WHERE toLower(d.name) CONTAINS toLower($q)
       RETURN d.id AS id, d.name AS name, d.categories AS categories,
              d.approved AS approved
       ORDER BY d.name
       LIMIT 20`,
      { q }
    );
    const drugs = result.records.map(r => ({
      id: r.get('id'),
      name: r.get('name'),
      categories: r.get('categories'),
      approved: r.get('approved'),
    }));
    res.json(drugs);
  } catch (err) {
    res.status(500).json({ error: err.message });
  } finally {
    await session.close();
  }
});

// GET /api/drugs/:id
router.get('/:id', async (req, res) => {
  const { id } = req.params;
  const session = openSession();
  try {
    const result = await session.run(
      `MATCH (d:Drug {id: $id})
       OPTIONAL MATCH (d)-[r:INTERACTS_WITH]-(other:Drug)
       OPTIONAL MATCH (d)-[:TARGETS]->(t:ProteinTarget)
       OPTIONAL MATCH (d)-[:PART_OF]->(p:Pathway)
       RETURN d,
         collect(DISTINCT {drug: other.name, drugId: other.id, severity: r.severity, description: r.description}) AS interactions,
         collect(DISTINCT {name: t.name, gene: t.gene}) AS targets,
         collect(DISTINCT p.name)[0..10] AS pathways`,
      { id }
    );
    if (result.records.length === 0) {
      return res.status(404).json({ error: 'Drug not found' });
    }
    const record = result.records[0];
    const drugNode = record.get('d').properties;
    const interactions = record.get('interactions').filter(i => i.drug !== null);
    const targets = record.get('targets').filter(t => t.name !== null);
    const pathways = record.get('pathways').filter(p => p !== null);

    res.json({ ...drugNode, interactions, targets, pathways });
  } catch (err) {
    res.status(500).json({ error: err.message });
  } finally {
    await session.close();
  }
});

// POST /api/drugs — Create Drug
router.post('/', async (req, res) => {
  const { id, name, description, categories, approved, drugType } = req.body;
  if (!id || !name) {
    return res.status(400).json({ error: 'id and name are required' });
  }
  const session = openSession();
  try {
    const result = await session.run(
      `CREATE (d:Drug {
        id: $id,
        name: $name,
        description: $description,
        categories: $categories,
        approved: $approved,
        drugType: $drugType
      })
      RETURN d`,
      { id, name, description: description || '', categories: categories || '', approved: approved !== false, drugType: drugType || '' }
    );
    const drug = result.records[0].get('d').properties;
    res.status(201).json(drug);
  } catch (err) {
    if (err.code === 'Neo.ClientError.Schema.ConstraintValidationFailed') {
      return res.status(409).json({ error: 'A drug with this ID already exists' });
    }
    res.status(500).json({ error: err.message });
  } finally {
    await session.close();
  }
});

// PUT /api/drugs/:id — Update Drug
router.put('/:id', async (req, res) => {
  const { id } = req.params;
  const { name, description, categories } = req.body;
  const session = openSession();
  try {
    const result = await session.run(
      `MATCH (d:Drug {id: $id})
       SET d.name = coalesce($name, d.name),
           d.description = coalesce($description, d.description),
           d.categories = coalesce($categories, d.categories)
       RETURN d`,
      { id, name: name || null, description: description || null, categories: categories || null }
    );
    if (result.records.length === 0) {
      return res.status(404).json({ error: 'Drug not found' });
    }
    res.json(result.records[0].get('d').properties);
  } catch (err) {
    res.status(500).json({ error: err.message });
  } finally {
    await session.close();
  }
});

// DELETE /api/drugs/:id
router.delete('/:id', async (req, res) => {
  const { id } = req.params;
  const session = openSession();
  try {
    await session.run(
      `MATCH (d:Drug {id: $id}) DETACH DELETE d`,
      { id }
    );
    res.json({ deleted: true });
  } catch (err) {
    res.status(500).json({ error: err.message });
  } finally {
    await session.close();
  }
});

export default router;
