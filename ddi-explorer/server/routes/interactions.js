import { Router } from 'express';
import driver, { session as openSession } from '../db.js';

const router = Router();
const VALID_SEVERITIES = ['Major', 'Moderate', 'Minor'];

// POST /api/interactions/check — Multi-Drug Checker
router.post('/check', async (req, res) => {
  const { drugIds } = req.body;
  if (!Array.isArray(drugIds) || drugIds.length < 2) {
    return res.status(400).json({ error: 'At least 2 drugIds required' });
  }
  const session = openSession();
  try {
    const result = await session.run(
      `MATCH (a:Drug)-[r:INTERACTS_WITH]-(b:Drug)
       WHERE a.id IN $drugIds
       AND b.id IN $drugIds
       AND id(a) < id(b)
       RETURN a.name AS drug1, a.id AS drug1Id,
              b.name AS drug2, b.id AS drug2Id,
              r.severity AS severity,
              r.description AS description
       ORDER BY
         CASE r.severity WHEN 'Major' THEN 1 WHEN 'Moderate' THEN 2 ELSE 3 END`,
      { drugIds }
    );
    const interactions = result.records.map(r => ({
      drug1: r.get('drug1'),
      drug1Id: r.get('drug1Id'),
      drug2: r.get('drug2'),
      drug2Id: r.get('drug2Id'),
      severity: r.get('severity'),
      description: r.get('description'),
    }));
    res.json(interactions);
  } catch (err) {
    res.status(500).json({ error: err.message });
  } finally {
    await session.close();
  }
});

// POST /api/interactions — Create Interaction
router.post('/', async (req, res) => {
  const { srcId, tgtId, severity, description } = req.body;
  if (!srcId || !tgtId || !severity) {
    return res.status(400).json({ error: 'srcId, tgtId, and severity are required' });
  }
  if (!VALID_SEVERITIES.includes(severity)) {
    return res.status(400).json({ error: 'severity must be Major, Moderate, or Minor' });
  }
  const session = openSession();
  try {
    const result = await session.run(
      `MATCH (a:Drug {id: $srcId}), (b:Drug {id: $tgtId})
       MERGE (a)-[r:INTERACTS_WITH]->(b)
       SET r.severity = $severity,
           r.description = $description
       RETURN r`,
      { srcId, tgtId, severity, description: description || '' }
    );
    if (result.records.length === 0) {
      return res.status(404).json({ error: 'One or both drugs not found' });
    }
    res.status(201).json(result.records[0].get('r').properties);
  } catch (err) {
    res.status(500).json({ error: err.message });
  } finally {
    await session.close();
  }
});

// PUT /api/interactions — Update Interaction
router.put('/', async (req, res) => {
  const { srcId, tgtId, severity, description } = req.body;
  if (!srcId || !tgtId || !severity) {
    return res.status(400).json({ error: 'srcId, tgtId, and severity are required' });
  }
  if (!VALID_SEVERITIES.includes(severity)) {
    return res.status(400).json({ error: 'severity must be Major, Moderate, or Minor' });
  }
  const session = openSession();
  try {
    const result = await session.run(
      `MATCH (a:Drug {id: $srcId})-[r:INTERACTS_WITH]->(b:Drug {id: $tgtId})
       SET r.severity = $severity,
           r.description = $description
       RETURN r`,
      { srcId, tgtId, severity, description: description || '' }
    );
    if (result.records.length === 0) {
      return res.status(404).json({ error: 'Interaction not found' });
    }
    res.json(result.records[0].get('r').properties);
  } catch (err) {
    res.status(500).json({ error: err.message });
  } finally {
    await session.close();
  }
});

// DELETE /api/interactions — Delete Interaction
router.delete('/', async (req, res) => {
  const { srcId, tgtId } = req.body;
  if (!srcId || !tgtId) {
    return res.status(400).json({ error: 'srcId and tgtId are required' });
  }
  const session = openSession();
  try {
    await session.run(
      `MATCH (a:Drug {id: $srcId})-[r:INTERACTS_WITH]->(b:Drug {id: $tgtId})
       DELETE r`,
      { srcId, tgtId }
    );
    res.json({ deleted: true });
  } catch (err) {
    res.status(500).json({ error: err.message });
  } finally {
    await session.close();
  }
});

export default router;
