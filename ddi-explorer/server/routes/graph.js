import { Router } from 'express';
import { session as openSession } from '../db.js';

const router = Router();

// GET /api/graph/:id?depth=1
router.get('/:id', async (req, res) => {
  const { id } = req.params;
  let depth = parseInt(req.query.depth) || 1;
  if (depth > 2) depth = 2;

  const session = openSession({ fetchSize: 200 });
  try {
    const query = depth === 1
      ? `MATCH path = (d:Drug {id: $id})-[:INTERACTS_WITH*1..1]-(neighbor:Drug)
         WITH d, neighbor,
              [r IN relationships(path) | {severity: r.severity}][0] AS rel
         RETURN d.id AS sourceId, d.name AS sourceName,
                neighbor.id AS targetId, neighbor.name AS targetName,
                rel.severity AS severity
         LIMIT 100`
      : `MATCH path = (d:Drug {id: $id})-[:INTERACTS_WITH*1..2]-(neighbor:Drug)
         WITH d, neighbor,
              [r IN relationships(path) | {severity: r.severity}][0] AS rel
         RETURN d.id AS sourceId, d.name AS sourceName,
                neighbor.id AS targetId, neighbor.name AS targetName,
                rel.severity AS severity
         LIMIT 100`;

    const result = await session.run(query, { id });
    const edges = result.records.map(r => ({
      sourceId: r.get('sourceId'),
      sourceName: r.get('sourceName'),
      targetId: r.get('targetId'),
      targetName: r.get('targetName'),
      severity: r.get('severity'),
    }));
    res.json(edges);
  } catch (err) {
    res.status(500).json({ error: err.message });
  } finally {
    await session.close();
  }
});

export default router;
