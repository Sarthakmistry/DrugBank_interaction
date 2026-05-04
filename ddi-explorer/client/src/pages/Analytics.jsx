import React, { useState, useEffect } from 'react';
import { useNavigate } from 'react-router-dom';
import {
  BarChart, Bar, XAxis, YAxis, Tooltip, ResponsiveContainer, Cell,
  PieChart, Pie, Legend,
} from 'recharts';
import {
  getAnalyticsSummary, getTopTargets, getTopInteractors, getSeverityBreakdown,
} from '../api.js';

function Spinner() {
  return (
    <div className="flex justify-center py-12">
      <div className="w-6 h-6 border-2 border-slate-300 border-t-blue-600 rounded-full animate-spin" />
    </div>
  );
}

function SummaryCard({ label, value, color }) {
  return (
    <div className={`rounded-xl p-5 border ${color}`}>
      <p className="text-sm font-medium text-slate-500 mb-1">{label}</p>
      <p className="text-3xl font-bold text-slate-900">{value?.toLocaleString() ?? '—'}</p>
    </div>
  );
}

const SEVERITY_COLORS = {
  Major: '#ef4444',
  Moderate: '#f97316',
  Minor: '#94a3b8',
};

export default function Analytics() {
  const [summary, setSummary] = useState(null);
  const [targets, setTargets] = useState([]);
  const [interactors, setInteractors] = useState([]);
  const [severity, setSeverity] = useState([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState('');
  const navigate = useNavigate();

  useEffect(() => {
    const load = async () => {
      setLoading(true);
      setError('');
      try {
        const [s, t, i, sv] = await Promise.all([
          getAnalyticsSummary(),
          getTopTargets(),
          getTopInteractors(),
          getSeverityBreakdown(),
        ]);
        setSummary(s);
        setTargets(t);
        setInteractors(i);
        setSeverity(sv);
      } catch (err) {
        setError(err.message);
      } finally {
        setLoading(false);
      }
    };
    load();
  }, []);

  if (loading) return <div className="p-6"><Spinner /></div>;

  if (error) return (
    <div className="p-6">
      <div className="p-4 bg-red-50 border border-red-200 text-red-700 rounded">{error}</div>
    </div>
  );

  const totalSeverity = severity.reduce((acc, d) => acc + d.count, 0);
  const severityWithPct = severity.map(d => ({
    ...d,
    pct: totalSeverity > 0 ? ((d.count / totalSeverity) * 100).toFixed(1) : 0,
    fill: SEVERITY_COLORS[d.severity] || '#cbd5e1',
  }));

  return (
    <div className="p-6 max-w-6xl mx-auto space-y-8">
      <h1 className="text-2xl font-bold text-slate-900">Analytics</h1>

      {/* Summary cards */}
      <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
        <SummaryCard label="Total Drugs" value={summary?.drugs} color="border-blue-200 bg-blue-50" />
        <SummaryCard label="Total Interactions" value={summary?.interactions} color="border-orange-200 bg-orange-50" />
        <SummaryCard label="Protein Targets" value={summary?.targets} color="border-purple-200 bg-purple-50" />
        <SummaryCard label="Pathways" value={summary?.pathways} color="border-green-200 bg-green-50" />
      </div>

      <div className="grid grid-cols-1 xl:grid-cols-2 gap-6">
        {/* Top Protein Targets */}
        <div className="bg-white border border-slate-200 rounded-xl p-5">
          <h2 className="text-base font-semibold text-slate-800 mb-4">Top 15 Most Targeted Proteins</h2>
          {targets.length === 0 ? (
            <p className="text-sm text-slate-400 text-center py-8">No data</p>
          ) : (
            <ResponsiveContainer width="100%" height={380}>
              <BarChart
                data={targets}
                layout="vertical"
                margin={{ top: 0, right: 20, left: 140, bottom: 0 }}
              >
                <XAxis type="number" tick={{ fontSize: 11 }} />
                <YAxis
                  type="category"
                  dataKey="gene"
                  tick={{ fontSize: 11 }}
                  width={130}
                  tickFormatter={(gene, i) => {
                    const t = targets[i];
                    return t ? `${gene || '?'} — ${(t.protein || '').substring(0, 14)}` : gene;
                  }}
                />
                <Tooltip
                  content={({ active, payload }) => {
                    if (!active || !payload?.length) return null;
                    const d = payload[0].payload;
                    return (
                      <div className="bg-white border border-slate-200 rounded p-2 text-xs shadow">
                        <p className="font-medium">{d.protein}</p>
                        <p className="text-slate-500">Gene: {d.gene}</p>
                        <p className="text-blue-600">{d.drugCount} drugs</p>
                      </div>
                    );
                  }}
                />
                <Bar dataKey="drugCount" fill="#3b82f6" radius={[0, 3, 3, 0]} />
              </BarChart>
            </ResponsiveContainer>
          )}
        </div>

        {/* Top Interactors */}
        <div className="bg-white border border-slate-200 rounded-xl p-5">
          <h2 className="text-base font-semibold text-slate-800 mb-4">Top 15 Most Interactive Drugs</h2>
          {interactors.length === 0 ? (
            <p className="text-sm text-slate-400 text-center py-8">No data</p>
          ) : (
            <ResponsiveContainer width="100%" height={380}>
              <BarChart
                data={interactors}
                layout="vertical"
                margin={{ top: 0, right: 20, left: 140, bottom: 0 }}
              >
                <XAxis type="number" tick={{ fontSize: 11 }} />
                <YAxis
                  type="category"
                  dataKey="drug"
                  tick={{ fontSize: 11 }}
                  width={130}
                  tickFormatter={v => v.length > 18 ? v.substring(0, 16) + '…' : v}
                />
                <Tooltip
                  content={({ active, payload }) => {
                    if (!active || !payload?.length) return null;
                    const d = payload[0].payload;
                    return (
                      <div className="bg-white border border-slate-200 rounded p-2 text-xs shadow">
                        <p className="font-medium">{d.drug}</p>
                        <p className="text-orange-600">{d.interactionCount} interactions</p>
                      </div>
                    );
                  }}
                />
                <Bar
                  dataKey="interactionCount"
                  radius={[0, 3, 3, 0]}
                  style={{ cursor: 'pointer' }}
                  onClick={(data) => navigate(`/?drug=${data.id}`)}
                >
                  {interactors.map((entry, i) => (
                    <Cell key={i} fill="#f97316" />
                  ))}
                </Bar>
              </BarChart>
            </ResponsiveContainer>
          )}
        </div>
      </div>

      {/* Severity breakdown */}
      <div className="bg-white border border-slate-200 rounded-xl p-5">
        <h2 className="text-base font-semibold text-slate-800 mb-4">Interaction Severity Breakdown</h2>
        {severity.length === 0 ? (
          <p className="text-sm text-slate-400 text-center py-8">No data</p>
        ) : (
          <div className="flex items-center gap-8">
            <ResponsiveContainer width={280} height={280}>
              <PieChart>
                <Pie
                  data={severityWithPct}
                  dataKey="count"
                  nameKey="severity"
                  cx="50%"
                  cy="50%"
                  innerRadius={70}
                  outerRadius={120}
                  paddingAngle={2}
                >
                  {severityWithPct.map((entry, i) => (
                    <Cell key={i} fill={entry.fill} />
                  ))}
                </Pie>
                <Tooltip
                  content={({ active, payload }) => {
                    if (!active || !payload?.length) return null;
                    const d = payload[0].payload;
                    return (
                      <div className="bg-white border border-slate-200 rounded p-2 text-xs shadow">
                        <p className="font-medium">{d.severity}</p>
                        <p>{d.count.toLocaleString()} ({d.pct}%)</p>
                      </div>
                    );
                  }}
                />
              </PieChart>
            </ResponsiveContainer>

            <div className="space-y-3">
              {severityWithPct.map(d => (
                <div key={d.severity} className="flex items-center gap-3">
                  <span className="w-3 h-3 rounded-full shrink-0" style={{ background: d.fill }} />
                  <div>
                    <p className="text-sm font-medium text-slate-800">{d.severity}</p>
                    <p className="text-xs text-slate-500">
                      {d.count.toLocaleString()} ({d.pct}%)
                    </p>
                  </div>
                </div>
              ))}
            </div>
          </div>
        )}
      </div>
    </div>
  );
}
