const BASE = import.meta.env.VITE_API_URL || '/api';

async function handleResponse(res) {
  const data = await res.json();
  if (!res.ok) throw new Error(data.error || 'Request failed');
  return data;
}

// Drugs
export const searchDrugs = (q) =>
  fetch(`${BASE}/drugs/search?q=${encodeURIComponent(q)}`).then(handleResponse);

export const getDrug = (id) =>
  fetch(`${BASE}/drugs/${id}`).then(handleResponse);

export const createDrug = (body) =>
  fetch(`${BASE}/drugs`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(body),
  }).then(handleResponse);

export const updateDrug = (id, body) =>
  fetch(`${BASE}/drugs/${id}`, {
    method: 'PUT',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(body),
  }).then(handleResponse);

export const deleteDrug = (id) =>
  fetch(`${BASE}/drugs/${id}`, { method: 'DELETE' }).then(handleResponse);

// Interactions
export const checkInteractions = (drugIds) =>
  fetch(`${BASE}/interactions/check`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ drugIds }),
  }).then(handleResponse);

export const createInteraction = (body) =>
  fetch(`${BASE}/interactions`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(body),
  }).then(handleResponse);

export const updateInteraction = (body) =>
  fetch(`${BASE}/interactions`, {
    method: 'PUT',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(body),
  }).then(handleResponse);

export const deleteInteraction = (body) =>
  fetch(`${BASE}/interactions`, {
    method: 'DELETE',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(body),
  }).then(handleResponse);

// Graph
export const getGraph = (id, depth = 1) =>
  fetch(`${BASE}/graph/${id}?depth=${depth}`).then(handleResponse);

// Analytics
export const getAnalyticsSummary = () =>
  fetch(`${BASE}/analytics/summary`).then(handleResponse);

export const getTopTargets = () =>
  fetch(`${BASE}/analytics/top-targets`).then(handleResponse);

export const getTopInteractors = () =>
  fetch(`${BASE}/analytics/top-interactors`).then(handleResponse);

export const getSeverityBreakdown = () =>
  fetch(`${BASE}/analytics/severity-breakdown`).then(handleResponse);
