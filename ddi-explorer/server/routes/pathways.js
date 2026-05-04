import { Router } from 'express';
import driver, { session as openSession } from '../db.js';

const router = Router();

// GET /api/pathways?drugId=<id>
router.get('/', async (req, res) => {
  const { drugId } = req.query;
  if (!drugId) return res.json([]);

  const session = openSession();
  try {
    const result = await session.run(
      `MATCH (d:Drug {id: $drugId})-[:PART_OF]->(p:Pathway)
       RETURN p.id AS id, p.name AS name, p.category AS category
       LIMIT 20`,
      { drugId }
    );
    const pathways = result.records.map(r => ({
      id: r.get('id'),
      name: r.get('name'),
      category: r.get('category'),
    }));
    res.json(pathways);
  } catch (err) {
    res.status(500).json({ error: err.message });
  } finally {
    await session.close();
  }
});

export default router;
