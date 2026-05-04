import { Router } from 'express';
import driver, { session as openSession } from '../db.js';

const router = Router();

// GET /api/targets?drugId=<id>
router.get('/', async (req, res) => {
  const { drugId } = req.query;
  const session = openSession();
  try {
    let query, params;
    if (drugId) {
      query = `MATCH (d:Drug {id: $drugId})-[:TARGETS]->(t:ProteinTarget)
               RETURN t.id AS id, t.name AS name, t.gene AS gene, t.uniprotId AS uniprotId`;
      params = { drugId };
    } else {
      query = `MATCH (t:ProteinTarget) RETURN t.id AS id, t.name AS name, t.gene AS gene, t.uniprotId AS uniprotId LIMIT 50`;
      params = {};
    }
    const result = await session.run(query, params);
    const targets = result.records.map(r => ({
      id: r.get('id'),
      name: r.get('name'),
      gene: r.get('gene'),
      uniprotId: r.get('uniprotId'),
    }));
    res.json(targets);
  } catch (err) {
    res.status(500).json({ error: err.message });
  } finally {
    await session.close();
  }
});

export default router;
