import { Router } from 'express';
import driver, { session as openSession } from '../db.js';

const router = Router();

// GET /api/analytics/summary
router.get('/summary', async (req, res) => {
  const session = openSession();
  try {
    const result = await session.run(
      `MATCH (d:Drug) WITH count(d) AS drugs
       MATCH ()-[r:INTERACTS_WITH]-() WITH drugs, count(r)/2 AS interactions
       MATCH (t:ProteinTarget) WITH drugs, interactions, count(t) AS targets
       MATCH (p:Pathway) RETURN drugs, interactions, targets, count(p) AS pathways`
    );
    const r = result.records[0];
    res.json({
      drugs: r.get('drugs').toNumber(),
      interactions: r.get('interactions').toNumber(),
      targets: r.get('targets').toNumber(),
      pathways: r.get('pathways').toNumber(),
    });
  } catch (err) {
    res.status(500).json({ error: err.message });
  } finally {
    await session.close();
  }
});

// GET /api/analytics/top-targets
router.get('/top-targets', async (req, res) => {
  const session = openSession();
  try {
    const result = await session.run(
      `MATCH (d:Drug)-[:TARGETS]->(t:ProteinTarget)
       RETURN t.name AS protein, t.gene AS gene, count(d) AS drugCount
       ORDER BY drugCount DESC
       LIMIT 15`
    );
    const data = result.records.map(r => ({
      protein: r.get('protein'),
      gene: r.get('gene'),
      drugCount: r.get('drugCount').toNumber(),
    }));
    res.json(data);
  } catch (err) {
    res.status(500).json({ error: err.message });
  } finally {
    await session.close();
  }
});

// GET /api/analytics/top-interactors
router.get('/top-interactors', async (req, res) => {
  const session = openSession();
  try {
    const result = await session.run(
      `MATCH (d:Drug)-[:INTERACTS_WITH]-()
       RETURN d.name AS drug, d.id AS id, count(*) AS interactionCount
       ORDER BY interactionCount DESC
       LIMIT 15`
    );
    const data = result.records.map(r => ({
      drug: r.get('drug'),
      id: r.get('id'),
      interactionCount: r.get('interactionCount').toNumber(),
    }));
    res.json(data);
  } catch (err) {
    res.status(500).json({ error: err.message });
  } finally {
    await session.close();
  }
});

// GET /api/analytics/severity-breakdown
router.get('/severity-breakdown', async (req, res) => {
  const session = openSession();
  try {
    const result = await session.run(
      `MATCH ()-[r:INTERACTS_WITH]-()
       RETURN r.severity AS severity, count(r) AS count`
    );
    const data = result.records.map(r => ({
      severity: r.get('severity'),
      count: r.get('count').toNumber(),
    }));
    res.json(data);
  } catch (err) {
    res.status(500).json({ error: err.message });
  } finally {
    await session.close();
  }
});

export default router;
