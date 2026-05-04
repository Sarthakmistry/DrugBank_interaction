import React, { useState, useEffect, useRef, useCallback } from 'react';
import { useNavigate } from 'react-router-dom';
import { searchDrugs, getDrug, deleteDrug, updateDrug } from '../api.js';

function SeverityBadge({ severity }) {
  const cls =
    severity === 'Major'
      ? 'bg-red-100 text-red-800'
      : severity === 'Moderate'
      ? 'bg-orange-100 text-orange-800'
      : 'bg-slate-100 text-slate-600';
  return (
    <span className={`inline-block px-2 py-0.5 rounded text-xs font-medium ${cls}`}>
      {severity}
    </span>
  );
}

function Spinner() {
  return (
    <div className="flex justify-center py-8">
      <div className="w-6 h-6 border-2 border-slate-300 border-t-blue-600 rounded-full animate-spin" />
    </div>
  );
}

function EditModal({ drug, onClose, onSaved }) {
  const [form, setForm] = useState({
    name: drug.name || '',
    description: drug.description || '',
    categories: drug.categories || '',
  });
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');

  const handleSubmit = async (e) => {
    e.preventDefault();
    setLoading(true);
    setError('');
    try {
      const updated = await updateDrug(drug.id, form);
      onSaved(updated);
    } catch (err) {
      setError(err.message);
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="fixed inset-0 bg-black/50 flex items-center justify-center z-50">
      <div className="bg-white rounded-lg shadow-xl w-full max-w-md p-6">
        <h2 className="text-lg font-semibold mb-4">Edit Drug — {drug.id}</h2>
        {error && <div className="mb-3 p-2 bg-red-50 text-red-700 text-sm rounded">{error}</div>}
        <form onSubmit={handleSubmit} className="space-y-3">
          <div>
            <label className="block text-sm font-medium text-slate-700 mb-1">Name</label>
            <input
              className="w-full border rounded px-3 py-2 text-sm"
              value={form.name}
              onChange={e => setForm(f => ({ ...f, name: e.target.value }))}
            />
          </div>
          <div>
            <label className="block text-sm font-medium text-slate-700 mb-1">Categories</label>
            <input
              className="w-full border rounded px-3 py-2 text-sm"
              value={form.categories}
              onChange={e => setForm(f => ({ ...f, categories: e.target.value }))}
            />
          </div>
          <div>
            <label className="block text-sm font-medium text-slate-700 mb-1">Description</label>
            <textarea
              className="w-full border rounded px-3 py-2 text-sm"
              rows={4}
              value={form.description}
              onChange={e => setForm(f => ({ ...f, description: e.target.value }))}
            />
          </div>
          <div className="flex gap-2 justify-end pt-2">
            <button type="button" onClick={onClose} className="px-4 py-2 text-sm rounded border hover:bg-slate-50">Cancel</button>
            <button
              type="submit"
              disabled={loading}
              className="px-4 py-2 text-sm rounded bg-blue-600 text-white hover:bg-blue-700 disabled:opacity-50"
            >
              {loading ? 'Saving…' : 'Save'}
            </button>
          </div>
        </form>
      </div>
    </div>
  );
}

export default function DrugSearch() {
  const [query, setQuery] = useState('');
  const [results, setResults] = useState([]);
  const [searchLoading, setSearchLoading] = useState(false);
  const [selectedDrug, setSelectedDrug] = useState(null);
  const [detailLoading, setDetailLoading] = useState(false);
  const [error, setError] = useState('');
  const [editOpen, setEditOpen] = useState(false);
  const [deleteConfirm, setDeleteConfirm] = useState(false);
  const [deleting, setDeleting] = useState(false);
  const debounceRef = useRef(null);
  const navigate = useNavigate();

  const doSearch = useCallback(async (q) => {
    if (q.length < 2) { setResults([]); return; }
    setSearchLoading(true);
    setError('');
    try {
      const data = await searchDrugs(q);
      setResults(data);
    } catch (err) {
      setError(err.message);
    } finally {
      setSearchLoading(false);
    }
  }, []);

  useEffect(() => {
    clearTimeout(debounceRef.current);
    debounceRef.current = setTimeout(() => doSearch(query), 300);
    return () => clearTimeout(debounceRef.current);
  }, [query, doSearch]);

  const selectDrug = async (drug) => {
    setDetailLoading(true);
    setError('');
    try {
      const detail = await getDrug(drug.id);
      setSelectedDrug(detail);
    } catch (err) {
      setError(err.message);
    } finally {
      setDetailLoading(false);
    }
  };

  const handleDelete = async () => {
    setDeleting(true);
    try {
      await deleteDrug(selectedDrug.id);
      setSelectedDrug(null);
      setDeleteConfirm(false);
      doSearch(query);
    } catch (err) {
      setError(err.message);
    } finally {
      setDeleting(false);
    }
  };

  return (
    <div className="p-6 max-w-6xl mx-auto">
      <h1 className="text-2xl font-bold text-slate-900 mb-6">Drug Search</h1>

      {error && (
        <div className="mb-4 p-3 bg-red-50 border border-red-200 text-red-700 rounded text-sm">{error}</div>
      )}

      <div className="flex gap-6">
        {/* Left panel */}
        <div className="w-80 shrink-0">
          <input
            type="text"
            placeholder="Search drugs (min 2 chars)…"
            className="w-full border border-slate-300 rounded-lg px-4 py-2.5 text-sm focus:outline-none focus:ring-2 focus:ring-blue-500"
            value={query}
            onChange={e => setQuery(e.target.value)}
          />

          <div className="mt-3 space-y-1">
            {searchLoading && <Spinner />}
            {!searchLoading && results.length === 0 && query.length >= 2 && (
              <p className="text-sm text-slate-500 py-4 text-center">No drugs found</p>
            )}
            {results.map(drug => (
              <button
                key={drug.id}
                onClick={() => selectDrug(drug)}
                className={`w-full text-left px-3 py-2.5 rounded-lg text-sm transition-colors ${
                  selectedDrug?.id === drug.id
                    ? 'bg-blue-50 border border-blue-200 text-blue-900'
                    : 'hover:bg-slate-100 text-slate-800'
                }`}
              >
                <div className="font-medium">{drug.name}</div>
                {drug.categories && (
                  <div className="text-xs text-slate-500 truncate mt-0.5">{drug.categories}</div>
                )}
              </button>
            ))}
          </div>
        </div>

        {/* Detail panel */}
        <div className="flex-1">
          {detailLoading && <Spinner />}

          {!detailLoading && selectedDrug && (
            <div className="bg-white border border-slate-200 rounded-xl shadow-sm overflow-hidden">
              {/* Header */}
              <div className="bg-[#0f172a] text-white px-6 py-4 flex items-start justify-between">
                <div>
                  <h2 className="text-xl font-semibold">{selectedDrug.name}</h2>
                  <p className="text-slate-400 text-sm mt-0.5">{selectedDrug.id}</p>
                </div>
                <div className="flex gap-2">
                  <button
                    onClick={() => setEditOpen(true)}
                    className="px-3 py-1.5 text-sm bg-slate-700 hover:bg-slate-600 rounded text-white"
                  >
                    Edit
                  </button>
                  <button
                    onClick={() => setDeleteConfirm(true)}
                    className="px-3 py-1.5 text-sm bg-red-700 hover:bg-red-600 rounded text-white"
                  >
                    Delete
                  </button>
                </div>
              </div>

              <div className="p-6 space-y-6">
                {/* Properties */}
                <section>
                  <h3 className="text-sm font-semibold text-slate-500 uppercase tracking-wide mb-3">Properties</h3>
                  <dl className="grid grid-cols-2 gap-x-4 gap-y-2 text-sm">
                    {[
                      ['Drug Type', selectedDrug.drugType],
                      ['Formula', selectedDrug.formula],
                      ['Mol Weight', selectedDrug.molWeight],
                      ['Half Life', selectedDrug.halfLife],
                      ['State', selectedDrug.state],
                      ['CAS Number', selectedDrug.casNumber],
                      ['Approved', selectedDrug.approved ? 'Yes' : 'No'],
                      ['Categories', selectedDrug.categories],
                    ].filter(([, v]) => v).map(([label, value]) => (
                      <div key={label}>
                        <dt className="text-slate-500">{label}</dt>
                        <dd className="text-slate-900 font-medium">{String(value)}</dd>
                      </div>
                    ))}
                  </dl>
                  {selectedDrug.indication && (
                    <div className="mt-3">
                      <dt className="text-slate-500 text-sm">Indication</dt>
                      <dd className="text-slate-800 text-sm mt-1">{selectedDrug.indication}</dd>
                    </div>
                  )}
                  {selectedDrug.mechanism && (
                    <div className="mt-3">
                      <dt className="text-slate-500 text-sm">Mechanism</dt>
                      <dd className="text-slate-800 text-sm mt-1">{selectedDrug.mechanism}</dd>
                    </div>
                  )}
                </section>

                {/* Interactions */}
                {selectedDrug.interactions?.length > 0 && (
                  <section>
                    <h3 className="text-sm font-semibold text-slate-500 uppercase tracking-wide mb-3">
                      Interactions ({selectedDrug.interactions.length})
                    </h3>
                    <div className="overflow-x-auto">
                      <table className="w-full text-sm">
                        <thead>
                          <tr className="border-b text-left text-slate-500">
                            <th className="pb-2 pr-4 font-medium">Drug</th>
                            <th className="pb-2 pr-4 font-medium">Severity</th>
                            <th className="pb-2 font-medium">Description</th>
                          </tr>
                        </thead>
                        <tbody className="divide-y divide-slate-100">
                          {selectedDrug.interactions.map((ix, i) => (
                            <tr key={i}>
                              <td className="py-2 pr-4">
                                <button
                                  className="text-blue-600 hover:underline text-left"
                                  onClick={() => selectDrug({ id: ix.drugId })}
                                >
                                  {ix.drug}
                                </button>
                              </td>
                              <td className="py-2 pr-4">
                                <SeverityBadge severity={ix.severity} />
                              </td>
                              <td className="py-2 text-slate-600 text-xs">{ix.description}</td>
                            </tr>
                          ))}
                        </tbody>
                      </table>
                    </div>
                  </section>
                )}

                {/* Targets */}
                {selectedDrug.targets?.length > 0 && (
                  <section>
                    <h3 className="text-sm font-semibold text-slate-500 uppercase tracking-wide mb-3">
                      Protein Targets ({selectedDrug.targets.length})
                    </h3>
                    <div className="flex flex-wrap gap-2">
                      {selectedDrug.targets.map((t, i) => (
                        <span key={i} className="inline-flex items-center gap-1 px-2 py-1 bg-blue-50 text-blue-800 rounded text-xs">
                          <span className="font-medium">{t.gene || '—'}</span>
                          <span className="text-blue-600">{t.name}</span>
                        </span>
                      ))}
                    </div>
                  </section>
                )}

                {/* Pathways */}
                {selectedDrug.pathways?.length > 0 && (
                  <section>
                    <h3 className="text-sm font-semibold text-slate-500 uppercase tracking-wide mb-3">
                      Pathways (first 10)
                    </h3>
                    <ul className="space-y-1">
                      {selectedDrug.pathways.map((p, i) => (
                        <li key={i} className="text-sm text-slate-700 flex items-center gap-1">
                          <span className="w-1.5 h-1.5 rounded-full bg-slate-400 shrink-0" />
                          {p}
                        </li>
                      ))}
                    </ul>
                  </section>
                )}
              </div>
            </div>
          )}

          {!detailLoading && !selectedDrug && (
            <div className="flex flex-col items-center justify-center h-64 text-slate-400">
              <p className="text-lg">Search for a drug to view details</p>
            </div>
          )}
        </div>
      </div>

      {/* Edit modal */}
      {editOpen && (
        <EditModal
          drug={selectedDrug}
          onClose={() => setEditOpen(false)}
          onSaved={(updated) => {
            setSelectedDrug(d => ({ ...d, ...updated }));
            setEditOpen(false);
          }}
        />
      )}

      {/* Delete confirm */}
      {deleteConfirm && (
        <div className="fixed inset-0 bg-black/50 flex items-center justify-center z-50">
          <div className="bg-white rounded-lg shadow-xl p-6 max-w-sm w-full">
            <h2 className="text-lg font-semibold mb-2">Delete Drug</h2>
            <p className="text-sm text-slate-600 mb-4">
              Are you sure you want to delete <strong>{selectedDrug.name}</strong>? This will also remove all its relationships.
            </p>
            <div className="flex gap-2 justify-end">
              <button onClick={() => setDeleteConfirm(false)} className="px-4 py-2 text-sm rounded border hover:bg-slate-50">Cancel</button>
              <button
                onClick={handleDelete}
                disabled={deleting}
                className="px-4 py-2 text-sm rounded bg-red-600 text-white hover:bg-red-700 disabled:opacity-50"
              >
                {deleting ? 'Deleting…' : 'Delete'}
              </button>
            </div>
          </div>
        </div>
      )}
    </div>
  );
}
