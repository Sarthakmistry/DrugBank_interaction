import React, { useState, useRef, useEffect, useCallback } from 'react';
import {
  searchDrugs, createDrug, updateDrug, deleteDrug,
  createInteraction, updateInteraction, deleteInteraction, getDrug,
} from '../api.js';

function Toast({ message, type, onClose }) {
  useEffect(() => {
    const t = setTimeout(onClose, 4000);
    return () => clearTimeout(t);
  }, [onClose]);
  return (
    <div className={`fixed bottom-6 right-6 z-50 px-4 py-3 rounded-lg shadow-lg text-sm text-white ${
      type === 'success' ? 'bg-green-600' : 'bg-red-600'
    }`}>
      {message}
    </div>
  );
}

function DrugAutocomplete({ label, onSelect, placeholder }) {
  const [q, setQ] = useState('');
  const [results, setResults] = useState([]);
  const [show, setShow] = useState(false);
  const debounceRef = useRef(null);

  useEffect(() => {
    clearTimeout(debounceRef.current);
    if (q.length < 2) { setResults([]); return; }
    debounceRef.current = setTimeout(async () => {
      try { setResults(await searchDrugs(q)); } catch { setResults([]); }
    }, 300);
  }, [q]);

  const select = (drug) => {
    setQ(drug.name);
    setShow(false);
    onSelect(drug);
  };

  return (
    <div className="relative">
      {label && <label className="block text-sm font-medium text-slate-700 mb-1">{label}</label>}
      <input
        type="text"
        placeholder={placeholder || 'Search drug…'}
        className="w-full border rounded px-3 py-2 text-sm focus:outline-none focus:ring-2 focus:ring-blue-500"
        value={q}
        onChange={e => { setQ(e.target.value); setShow(true); }}
        onFocus={() => setShow(true)}
        onBlur={() => setTimeout(() => setShow(false), 150)}
      />
      {show && results.length > 0 && (
        <div className="absolute top-full left-0 right-0 bg-white border border-slate-200 rounded shadow-lg mt-1 z-20 max-h-48 overflow-y-auto">
          {results.map(d => (
            <button
              key={d.id}
              onMouseDown={() => select(d)}
              className="w-full text-left px-3 py-2 text-sm hover:bg-slate-50"
            >
              {d.name} <span className="text-slate-400 text-xs">{d.id}</span>
            </button>
          ))}
        </div>
      )}
    </div>
  );
}

// --- Drugs Tab ---
function DrugsTab({ onToast }) {
  const [mode, setMode] = useState('create'); // create | edit | delete
  const [form, setForm] = useState({ id: '', name: '', description: '', categories: '', drugType: '', approved: true });
  const [idError, setIdError] = useState('');
  const [loading, setLoading] = useState(false);
  const [deleteConfirm, setDeleteConfirm] = useState(false);
  const [selectedDrug, setSelectedDrug] = useState(null);

  const validateId = (id) => /^DB\d{5}$/.test(id);

  const handleCreate = async (e) => {
    e.preventDefault();
    if (!validateId(form.id)) { setIdError('ID must match format DB00000'); return; }
    setIdError('');
    setLoading(true);
    try {
      await createDrug(form);
      onToast('Drug created successfully', 'success');
      setForm({ id: '', name: '', description: '', categories: '', drugType: '', approved: true });
    } catch (err) {
      onToast(err.message, 'error');
    } finally {
      setLoading(false);
    }
  };

  const loadForEdit = async (drug) => {
    setSelectedDrug(drug);
    try {
      const detail = await getDrug(drug.id);
      setForm({ name: detail.name || '', description: detail.description || '', categories: detail.categories || '' });
    } catch (err) {
      onToast(err.message, 'error');
    }
  };

  const handleEdit = async (e) => {
    e.preventDefault();
    if (!selectedDrug) return;
    setLoading(true);
    try {
      await updateDrug(selectedDrug.id, form);
      onToast('Drug updated successfully', 'success');
    } catch (err) {
      onToast(err.message, 'error');
    } finally {
      setLoading(false);
    }
  };

  const handleDelete = async () => {
    if (!selectedDrug) return;
    setLoading(true);
    try {
      await deleteDrug(selectedDrug.id);
      onToast('Drug deleted', 'success');
      setSelectedDrug(null);
      setDeleteConfirm(false);
    } catch (err) {
      onToast(err.message, 'error');
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="space-y-4">
      {/* Mode tabs */}
      <div className="flex gap-1 border-b border-slate-200 pb-0">
        {['create', 'edit', 'delete'].map(m => (
          <button
            key={m}
            onClick={() => setMode(m)}
            className={`px-4 py-2 text-sm capitalize transition-colors border-b-2 -mb-px ${
              mode === m ? 'border-blue-600 text-blue-600' : 'border-transparent text-slate-500 hover:text-slate-800'
            }`}
          >
            {m}
          </button>
        ))}
      </div>

      {mode === 'create' && (
        <form onSubmit={handleCreate} className="space-y-3 max-w-md">
          <div>
            <label className="block text-sm font-medium text-slate-700 mb-1">Drug ID</label>
            <input
              className={`w-full border rounded px-3 py-2 text-sm ${idError ? 'border-red-400' : ''}`}
              placeholder="DB00001"
              value={form.id}
              onChange={e => { setForm(f => ({ ...f, id: e.target.value })); setIdError(''); }}
            />
            {idError && <p className="text-red-600 text-xs mt-1">{idError}</p>}
          </div>
          <div>
            <label className="block text-sm font-medium text-slate-700 mb-1">Name *</label>
            <input
              className="w-full border rounded px-3 py-2 text-sm"
              value={form.name}
              onChange={e => setForm(f => ({ ...f, name: e.target.value }))}
              required
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
            <label className="block text-sm font-medium text-slate-700 mb-1">Drug Type</label>
            <input
              className="w-full border rounded px-3 py-2 text-sm"
              value={form.drugType}
              onChange={e => setForm(f => ({ ...f, drugType: e.target.value }))}
            />
          </div>
          <div>
            <label className="block text-sm font-medium text-slate-700 mb-1">Description</label>
            <textarea
              className="w-full border rounded px-3 py-2 text-sm"
              rows={3}
              value={form.description}
              onChange={e => setForm(f => ({ ...f, description: e.target.value }))}
            />
          </div>
          <label className="flex items-center gap-2 text-sm">
            <input type="checkbox" checked={form.approved} onChange={e => setForm(f => ({ ...f, approved: e.target.checked }))} />
            Approved
          </label>
          <button
            type="submit"
            disabled={loading}
            className="px-5 py-2 bg-blue-600 text-white text-sm rounded hover:bg-blue-700 disabled:opacity-50"
          >
            {loading ? 'Creating…' : 'Create Drug'}
          </button>
        </form>
      )}

      {mode === 'edit' && (
        <div className="space-y-4 max-w-md">
          <DrugAutocomplete label="Search drug to edit" onSelect={loadForEdit} />
          {selectedDrug && (
            <form onSubmit={handleEdit} className="space-y-3">
              <p className="text-xs text-slate-500">Editing: <strong>{selectedDrug.id}</strong></p>
              {[['Name', 'name'], ['Categories', 'categories']].map(([lbl, key]) => (
                <div key={key}>
                  <label className="block text-sm font-medium text-slate-700 mb-1">{lbl}</label>
                  <input
                    className="w-full border rounded px-3 py-2 text-sm"
                    value={form[key] || ''}
                    onChange={e => setForm(f => ({ ...f, [key]: e.target.value }))}
                  />
                </div>
              ))}
              <div>
                <label className="block text-sm font-medium text-slate-700 mb-1">Description</label>
                <textarea
                  className="w-full border rounded px-3 py-2 text-sm"
                  rows={3}
                  value={form.description || ''}
                  onChange={e => setForm(f => ({ ...f, description: e.target.value }))}
                />
              </div>
              <button
                type="submit"
                disabled={loading}
                className="px-5 py-2 bg-blue-600 text-white text-sm rounded hover:bg-blue-700 disabled:opacity-50"
              >
                {loading ? 'Saving…' : 'Save Changes'}
              </button>
            </form>
          )}
        </div>
      )}

      {mode === 'delete' && (
        <div className="space-y-4 max-w-md">
          <DrugAutocomplete label="Search drug to delete" onSelect={d => { setSelectedDrug(d); setDeleteConfirm(false); }} />
          {selectedDrug && (
            <div>
              <p className="text-sm text-slate-700 mb-3">
                Selected: <strong>{selectedDrug.name}</strong> ({selectedDrug.id})
              </p>
              {!deleteConfirm ? (
                <button
                  onClick={() => setDeleteConfirm(true)}
                  className="px-4 py-2 bg-red-600 text-white text-sm rounded hover:bg-red-700"
                >
                  Delete Drug
                </button>
              ) : (
                <div className="bg-red-50 border border-red-200 rounded p-4">
                  <p className="text-sm text-red-800 mb-3">
                    This will permanently delete <strong>{selectedDrug.name}</strong> and all its relationships.
                  </p>
                  <div className="flex gap-2">
                    <button
                      onClick={() => setDeleteConfirm(false)}
                      className="px-4 py-2 text-sm border rounded hover:bg-slate-50"
                    >
                      Cancel
                    </button>
                    <button
                      onClick={handleDelete}
                      disabled={loading}
                      className="px-4 py-2 text-sm bg-red-600 text-white rounded hover:bg-red-700 disabled:opacity-50"
                    >
                      {loading ? 'Deleting…' : 'Confirm Delete'}
                    </button>
                  </div>
                </div>
              )}
            </div>
          )}
        </div>
      )}
    </div>
  );
}

// --- Interactions Tab ---
function InteractionsTab({ onToast }) {
  const [mode, setMode] = useState('create');
  const [src, setSrc] = useState(null);
  const [tgt, setTgt] = useState(null);
  const [severity, setSeverity] = useState('Major');
  const [description, setDescription] = useState('');
  const [loading, setLoading] = useState(false);
  const [deleteConfirm, setDeleteConfirm] = useState(false);

  const handleCreate = async (e) => {
    e.preventDefault();
    if (!src || !tgt) { onToast('Select both source and target drugs', 'error'); return; }
    setLoading(true);
    try {
      await createInteraction({ srcId: src.id, tgtId: tgt.id, severity, description });
      onToast('Interaction created', 'success');
      setSrc(null); setTgt(null); setDescription('');
    } catch (err) {
      onToast(err.message, 'error');
    } finally {
      setLoading(false);
    }
  };

  const handleUpdate = async (e) => {
    e.preventDefault();
    if (!src || !tgt) { onToast('Select both source and target drugs', 'error'); return; }
    setLoading(true);
    try {
      await updateInteraction({ srcId: src.id, tgtId: tgt.id, severity, description });
      onToast('Interaction updated', 'success');
    } catch (err) {
      onToast(err.message, 'error');
    } finally {
      setLoading(false);
    }
  };

  const handleDelete = async () => {
    if (!src || !tgt) return;
    setLoading(true);
    try {
      await deleteInteraction({ srcId: src.id, tgtId: tgt.id });
      onToast('Interaction deleted', 'success');
      setSrc(null); setTgt(null); setDeleteConfirm(false);
    } catch (err) {
      onToast(err.message, 'error');
    } finally {
      setLoading(false);
    }
  };

  const drugSelectors = (
    <div className="grid grid-cols-2 gap-4">
      <DrugAutocomplete label="Source Drug" onSelect={setSrc} placeholder="Search source…" />
      <DrugAutocomplete label="Target Drug" onSelect={setTgt} placeholder="Search target…" />
    </div>
  );

  const severityField = (
    <div>
      <label className="block text-sm font-medium text-slate-700 mb-1">Severity</label>
      <select
        className="w-full border rounded px-3 py-2 text-sm"
        value={severity}
        onChange={e => setSeverity(e.target.value)}
      >
        <option>Major</option>
        <option>Moderate</option>
        <option>Minor</option>
      </select>
    </div>
  );

  const descriptionField = (
    <div>
      <label className="block text-sm font-medium text-slate-700 mb-1">Description</label>
      <textarea
        className="w-full border rounded px-3 py-2 text-sm"
        rows={3}
        value={description}
        onChange={e => setDescription(e.target.value)}
      />
    </div>
  );

  return (
    <div className="space-y-4">
      <div className="flex gap-1 border-b border-slate-200 pb-0">
        {['create', 'edit', 'delete'].map(m => (
          <button
            key={m}
            onClick={() => setMode(m)}
            className={`px-4 py-2 text-sm capitalize transition-colors border-b-2 -mb-px ${
              mode === m ? 'border-blue-600 text-blue-600' : 'border-transparent text-slate-500 hover:text-slate-800'
            }`}
          >
            {m}
          </button>
        ))}
      </div>

      {mode === 'create' && (
        <form onSubmit={handleCreate} className="space-y-4 max-w-lg">
          {drugSelectors}
          {severityField}
          {descriptionField}
          <button
            type="submit"
            disabled={loading}
            className="px-5 py-2 bg-blue-600 text-white text-sm rounded hover:bg-blue-700 disabled:opacity-50"
          >
            {loading ? 'Creating…' : 'Create Interaction'}
          </button>
        </form>
      )}

      {mode === 'edit' && (
        <form onSubmit={handleUpdate} className="space-y-4 max-w-lg">
          {drugSelectors}
          {severityField}
          {descriptionField}
          <button
            type="submit"
            disabled={loading}
            className="px-5 py-2 bg-blue-600 text-white text-sm rounded hover:bg-blue-700 disabled:opacity-50"
          >
            {loading ? 'Saving…' : 'Update Interaction'}
          </button>
        </form>
      )}

      {mode === 'delete' && (
        <div className="space-y-4 max-w-lg">
          {drugSelectors}
          {src && tgt && !deleteConfirm && (
            <button
              onClick={() => setDeleteConfirm(true)}
              className="px-4 py-2 bg-red-600 text-white text-sm rounded hover:bg-red-700"
            >
              Delete Interaction
            </button>
          )}
          {deleteConfirm && (
            <div className="bg-red-50 border border-red-200 rounded p-4">
              <p className="text-sm text-red-800 mb-3">
                Delete interaction between <strong>{src?.name}</strong> and <strong>{tgt?.name}</strong>?
              </p>
              <div className="flex gap-2">
                <button onClick={() => setDeleteConfirm(false)} className="px-4 py-2 text-sm border rounded hover:bg-slate-50">Cancel</button>
                <button
                  onClick={handleDelete}
                  disabled={loading}
                  className="px-4 py-2 text-sm bg-red-600 text-white rounded hover:bg-red-700 disabled:opacity-50"
                >
                  {loading ? 'Deleting…' : 'Confirm Delete'}
                </button>
              </div>
            </div>
          )}
        </div>
      )}
    </div>
  );
}

export default function CRUDManager() {
  const [tab, setTab] = useState('drugs');
  const [toast, setToast] = useState(null);

  const showToast = (message, type) => setToast({ message, type });

  return (
    <div className="p-6 max-w-3xl mx-auto">
      <h1 className="text-2xl font-bold text-slate-900 mb-6">CRUD Manager</h1>

      {/* Main tabs */}
      <div className="flex gap-0 mb-6 border border-slate-200 rounded-lg overflow-hidden w-fit">
        {['drugs', 'interactions'].map(t => (
          <button
            key={t}
            onClick={() => setTab(t)}
            className={`px-6 py-2.5 text-sm font-medium capitalize transition-colors ${
              tab === t ? 'bg-[#0f172a] text-white' : 'text-slate-600 hover:bg-slate-50'
            }`}
          >
            {t}
          </button>
        ))}
      </div>

      <div className="bg-white border border-slate-200 rounded-xl p-6">
        {tab === 'drugs' ? (
          <DrugsTab onToast={showToast} />
        ) : (
          <InteractionsTab onToast={showToast} />
        )}
      </div>

      {toast && (
        <Toast
          message={toast.message}
          type={toast.type}
          onClose={() => setToast(null)}
        />
      )}
    </div>
  );
}
