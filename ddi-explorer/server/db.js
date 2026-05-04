import neo4j from 'neo4j-driver';
import dotenv from 'dotenv';
dotenv.config({ path: '../.env' });

const driver = neo4j.driver(
  process.env.NEO4J_URI,
  neo4j.auth.basic(process.env.NEO4J_USERNAME, process.env.NEO4J_PASSWORD)
);

// Helper: open a session targeting the correct database.
export function session(config = {}) {
  return driver.session({
    database: process.env.NEO4J_DATABASE || 'neo4j',
    ...config,
  });
}

export default driver;
