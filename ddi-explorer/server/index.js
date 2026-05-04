import express from 'express';
import cors from 'cors';
import dotenv from 'dotenv';
dotenv.config({ path: '../.env' });

import drugsRouter from './routes/drugs.js';
import interactionsRouter from './routes/interactions.js';
import targetsRouter from './routes/targets.js';
import pathwaysRouter from './routes/pathways.js';
import analyticsRouter from './routes/analytics.js';

const app = express();
const PORT = process.env.PORT || 3001;

app.use(cors());
app.use(express.json());

app.use('/api/drugs', drugsRouter);
app.use('/api/interactions', interactionsRouter);
app.use('/api/targets', targetsRouter);
app.use('/api/pathways', pathwaysRouter);
app.use('/api/analytics', analyticsRouter);

// Graph endpoint
import graphRouter from './routes/graph.js';
app.use('/api/graph', graphRouter);

// Error handling middleware
app.use((err, req, res, next) => {
  console.error(err.stack);
  res.status(500).json({ error: err.message || 'Internal server error' });
});

app.listen(PORT, () => {
  console.log(`DDI Explorer API running on http://localhost:${PORT}`);
});
