import React, { useState, useEffect, useRef, useCallback } from 'react';
import { searchDrugs, checkInteractions } from '../api.js';

function SeverityBadge({ severity }) {
  const cls =
    severity === 'Major'
      ? 'bg-red-100 text-red-800'
      : severity === 'Moderate'
      ? 'bg-orange-100 text-orange-800'
      : 'bg-slate-100 text-slate-600';
  return (
    <span className={`inline-block px-2 py-0.5 rounded-full text-xs font-medium ${cls}`}>
      {severity}
    </span>
  );
}

export default function InteractionChecker() {
  const [query, setQuery] = useState('');
  const [suggestions, setSuggestions] = useState([]);
  const [regimen, setRegimen] = useState([]);
  const [interactions, setInteractions] = useState(null);
  const [filter, setFilter] = useState('all');
  const [loading, setLoading] = useState(false);
  const [checking, setChecking] = useState(false);
  const [error, setError] = useState('');
  const [showSuggestions, setShowSuggestions] = useState(false);
  const debounceRef = useRef(null);

  const doSearch = useCallback(async (q) => {
    if (q.length < 2) { setSuggestions([]); return; }
    setLoading(true);
    try {
      const data = await searchDrugs(q);
      setSuggestions(data.filter(d => !regimen.find(r => r.id === d.id)));
    } catch {
      setSuggestions([]);
    } finally {
      setLoading(false);
    }
  }, [regimen]);

  useEffect(() => {
    clearTimeout(debounceRef.current);
    debounceRef.current = setTimeout(() => doSearch(query), 300);
    return () => clearTimeout(debounceRef.current);
  }, [query, doSearch]);

  const addDrug = (drug) => {
    if (regimen.length >= 10) return;
    if (regimen.find(d => d.id === drug.id)) return;
    setRegimen(r => [...r, drug]);
    setQuery('');
    setSuggestions([]);
    setShowSuggestions(false);
    setInteractions(null);
  };

  const removeDrug = (id) => {
    setRegimen(r => r.filter(d => d.id !== id));
    setInteractions(null);
  };

  const handleCheck = async () => {
    setChecking(true);
    setError('');
    try {
      const data = await checkInteractions(regimen.map(d => d.id));
      setInteractions(data);
    } catch (err) {
      setError(err.message);
    } finally {
      setChecking(false);
    }
  };

  const filtered = interactions?.filter(ix => {
    if (filter === 'major') return ix.severity === 'Major';
    if (filter === 'major-moderate') return ['Major', 'Moderate'].includes(ix.severity);
    return true;
  }) ?? [];

  const majorCount = interactions?.filter(i => i.severity === 'Major').length ?? 0;

  return (
    <div className="p-6 max-w-4xl mx-auto">
      <h1 className="text-2xl font-bold text-slate-900 mb-6">Interaction Checker</h1>

      {error && (
        <div className="mb-4 p-3 bg-red-50 border border-red-200 text-red-700 rounded text-sm">{error}</div>
      )}

      {/* Drug search input */}
      <div className="relative mb-4">
        <input
          type="text"
          placeholder={regimen.length >= 10 ? 'Maximum 10 drugs reached' : 'Add drug to regimen…'}
          className="w-full border border-slate-300 rounded-lg px-4 py-2.5 text-sm focus:outline-none focus:ring-2 focus:ring-blue-500"
          value={query}
          disabled={regimen.length >= 10}
          onChange={e => { setQuery(e.target.value); setShowSuggestions(true); }}
          onFocus={() => setShowSuggestions(true)}
          onBlur={() => setTimeout(() => setShowSuggestions(false), 150)}
        />
        {showSuggestions && suggestions.length > 0 && (
          <div className="absolute top-full left-0 right-0 bg-white border border-slate-200 rounded-lg shadow-lg mt-1 z-10 max-h-60 overflow-y-auto">
            {suggestions.map(drug => (
              <button
                key={drug.id}
                onMouseDown={() => addDrug(drug)}
                className="w-full text-left px-4 py-2 hover:bg-slate-50 text-sm"
              >
                <span className="font-medium">{drug.name}</span>
                <span className="text-slate-400 ml-2 text-xs">{drug.id}</span>
              </button>
            ))}
          </div>
        )}
      </div>

      {/* Regimen list */}
      {regimen.length > 0 && (
        <div className="flex flex-wrap gap-2 mb-6">
          {regimen.map(drug => (
            <span
              key={drug.id}
              className="inline-flex items-center gap-1.5 px-3 py-1.5 bg-[#0f172a] text-white rounded-full text-sm"
            >
              {drug.name}
              <button
                onClick={() => removeDrug(drug.id)}
                className="hover:text-red-300 ml-1 text-slate-400 text-xs leading-none"
              >
                ×
              </button>
            </span>
          ))}
        </div>
      )}

      {/* Check button */}
      <button
        onClick={handleCheck}
        disabled={regimen.length < 2 || checking}
        className="px-6 py-2.5 bg-blue-600 text-white rounded-lg text-sm font-medium hover:bg-blue-700 disabled:opacity-40 disabled:cursor-not-allowed"
      >
        {checking ? 'Checking…' : 'Check Interactions'}
      </button>

      {/* Results */}
      {interactions !== null && (
        <div className="mt-6">
          {interactions.length === 0 ? (
            <div className="bg-green-50 border border-green-200 text-green-800 rounded-lg p-4 text-sm font-medium">
              No known interactions found among {regimen.length} drugs.
            </div>
          ) : (
            <>
              {/* Summary + filter */}
              <div className="flex items-center justify-between mb-4">
                <p className="text-sm text-slate-700">
                  <strong>{interactions.length}</strong> interaction{interactions.length !== 1 ? 's' : ''} found among{' '}
                  <strong>{regimen.length}</strong> drugs
                  {majorCount > 0 && (
                    <span className="ml-1 text-red-700 font-medium">({majorCount} major)</span>
                  )}
                </p>
                <div className="flex gap-1 text-xs">
                  {['all', 'major', 'major-moderate'].map(f => (
                    <button
                      key={f}
                      onClick={() => setFilter(f)}
                      className={`px-2.5 py-1 rounded border transition-colors ${
                        filter === f ? 'bg-slate-800 text-white border-slate-800' : 'text-slate-600 hover:bg-slate-50'
                      }`}
                    >
                      {f === 'all' ? 'All' : f === 'major' ? 'Major only' : 'Major + Moderate'}
                    </button>
                  ))}
                </div>
              </div>

              {filtered.length === 0 ? (
                <p className="text-sm text-slate-500 text-center py-4">No interactions match this filter.</p>
              ) : (
                <div className="overflow-x-auto rounded-lg border border-slate-200">
                  <table className="w-full text-sm">
                    <thead className="bg-slate-50 border-b border-slate-200">
                      <tr>
                        <th className="text-left px-4 py-3 text-slate-600 font-medium">Drug 1</th>
                        <th className="text-left px-4 py-3 text-slate-600 font-medium">Drug 2</th>
                        <th className="text-left px-4 py-3 text-slate-600 font-medium">Severity</th>
                        <th className="text-left px-4 py-3 text-slate-600 font-medium">Description</th>
                      </tr>
                    </thead>
                    <tbody className="divide-y divide-slate-100">
                      {filtered.map((ix, i) => (
                        <tr key={i} className="hover:bg-slate-50">
                          <td className="px-4 py-3 font-medium">{ix.drug1}</td>
                          <td className="px-4 py-3 font-medium">{ix.drug2}</td>
                          <td className="px-4 py-3">
                            <SeverityBadge severity={ix.severity} />
                          </td>
                          <td className="px-4 py-3 text-slate-600 text-xs max-w-xs">{ix.description}</td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              )}
            </>
          )}
        </div>
      )}
    </div>
  );
}
