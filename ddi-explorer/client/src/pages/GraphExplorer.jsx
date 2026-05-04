import React, { useState, useEffect, useRef, useCallback } from 'react';
import * as d3 from 'd3';
import { searchDrugs, getGraph } from '../api.js';

const SEVERITY_COLOR = {
  Major: '#ef4444',
  Moderate: '#f97316',
  Minor: '#94a3b8',
};

export default function GraphExplorer() {
  const svgRef = useRef(null);
  const simulationRef = useRef(null);

  const [query, setQuery] = useState('');
  const [suggestions, setSuggestions] = useState([]);
  const [showSugg, setShowSugg] = useState(false);
  const [rootDrug, setRootDrug] = useState(null);
  const [depth, setDepth] = useState(1);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');
  const [tooltip, setTooltip] = useState(null);
  const debounceRef = useRef(null);

  // Search suggestions
  const doSearch = useCallback(async (q) => {
    if (q.length < 2) { setSuggestions([]); return; }
    try {
      const data = await searchDrugs(q);
      setSuggestions(data);
    } catch { setSuggestions([]); }
  }, []);

  useEffect(() => {
    clearTimeout(debounceRef.current);
    debounceRef.current = setTimeout(() => doSearch(query), 300);
    return () => clearTimeout(debounceRef.current);
  }, [query, doSearch]);

  // Build graph
  const buildGraph = useCallback(async (drug) => {
    setLoading(true);
    setError('');
    try {
      const edges = await getGraph(drug.id, depth);

      if (simulationRef.current) simulationRef.current.stop();
      const svg = d3.select(svgRef.current);
      svg.selectAll('*').remove();

      const width = svgRef.current.clientWidth || 900;
      const height = svgRef.current.clientHeight || 600;

      // Build nodes/links from edges
      const nodeMap = new Map();
      const links = [];

      edges.forEach(e => {
        if (!nodeMap.has(e.sourceId)) {
          nodeMap.set(e.sourceId, { id: e.sourceId, name: e.sourceName, isRoot: e.sourceId === drug.id });
        }
        if (!nodeMap.has(e.targetId)) {
          nodeMap.set(e.targetId, { id: e.targetId, name: e.targetName, isRoot: false });
        }
        links.push({ source: e.sourceId, target: e.targetId, severity: e.severity });
      });

      // Compute worst severity per node
      const worstSeverity = new Map();
      links.forEach(l => {
        const update = (id) => {
          const cur = worstSeverity.get(id);
          if (!cur || (l.severity === 'Major') || (l.severity === 'Moderate' && cur !== 'Major')) {
            worstSeverity.set(id, l.severity);
          }
        };
        update(typeof l.source === 'object' ? l.source.id : l.source);
        update(typeof l.target === 'object' ? l.target.id : l.target);
      });

      // Count connections per node
      const connectionCount = new Map();
      links.forEach(l => {
        const sid = typeof l.source === 'object' ? l.source.id : l.source;
        const tid = typeof l.target === 'object' ? l.target.id : l.target;
        connectionCount.set(sid, (connectionCount.get(sid) || 0) + 1);
        connectionCount.set(tid, (connectionCount.get(tid) || 0) + 1);
      });

      const nodes = Array.from(nodeMap.values());

      const zoom = d3.zoom().scaleExtent([0.2, 4]).on('zoom', (event) => {
        g.attr('transform', event.transform);
      });

      svg.call(zoom);
      const g = svg.append('g');

      // Defs for arrow markers
      svg.append('defs').selectAll('marker')
        .data(['Major', 'Moderate', 'Minor'])
        .join('marker')
        .attr('id', d => `arrow-${d}`)
        .attr('viewBox', '0 -5 10 10')
        .attr('refX', 20)
        .attr('refY', 0)
        .attr('markerWidth', 6)
        .attr('markerHeight', 6)
        .attr('orient', 'auto')
        .append('path')
        .attr('d', 'M0,-5L10,0L0,5')
        .attr('fill', d => SEVERITY_COLOR[d] || '#94a3b8');

      const simulation = d3.forceSimulation(nodes)
        .force('link', d3.forceLink(links).id(d => d.id).distance(80))
        .force('charge', d3.forceManyBody().strength(-300))
        .force('center', d3.forceCenter(width / 2, height / 2))
        .force('collision', d3.forceCollide(30));

      simulationRef.current = simulation;

      const link = g.append('g').selectAll('line')
        .data(links)
        .join('line')
        .attr('stroke', d => SEVERITY_COLOR[d.severity] || '#94a3b8')
        .attr('stroke-width', d => d.severity === 'Major' ? 3 : d.severity === 'Moderate' ? 2 : 1)
        .attr('stroke-opacity', 0.7)
        .attr('marker-end', d => `url(#arrow-${d.severity || 'Minor'})`)
        .style('cursor', 'pointer')
        .on('mouseover', (event, d) => {
          const desc = d.description ? d.description.substring(0, 120) + (d.description.length > 120 ? '…' : '') : 'No description';
          setTooltip({ x: event.pageX, y: event.pageY, content: `${d.severity || 'Unknown'}: ${desc}` });
        })
        .on('mouseout', () => setTooltip(null));

      const node = g.append('g').selectAll('circle')
        .data(nodes)
        .join('circle')
        .attr('r', d => d.isRoot ? 18 : 10)
        .attr('fill', d => {
          if (d.isRoot) return '#0f172a';
          const sev = worstSeverity.get(d.id);
          return SEVERITY_COLOR[sev] || '#64748b';
        })
        .attr('stroke', '#fff')
        .attr('stroke-width', 2)
        .style('cursor', 'pointer')
        .on('mouseover', (event, d) => {
          const cnt = connectionCount.get(d.id) || 0;
          setTooltip({ x: event.pageX, y: event.pageY, content: `${d.name} — ${cnt} connection${cnt !== 1 ? 's' : ''}` });
        })
        .on('mouseout', () => setTooltip(null))
        .on('click', (event, d) => {
          setRootDrug({ id: d.id, name: d.name });
        })
        .call(
          d3.drag()
            .on('start', (event, d) => {
              if (!event.active) simulation.alphaTarget(0.3).restart();
              d.fx = d.x; d.fy = d.y;
            })
            .on('drag', (event, d) => { d.fx = event.x; d.fy = event.y; })
            .on('end', (event, d) => {
              if (!event.active) simulation.alphaTarget(0);
              d.fx = null; d.fy = null;
            })
        );

      // Entry animation
      node.attr('r', 0).transition().duration(400).attr('r', d => d.isRoot ? 18 : 10);

      const label = g.append('g').selectAll('text')
        .data(nodes)
        .join('text')
        .attr('font-size', d => d.isRoot ? 13 : 10)
        .attr('font-weight', d => d.isRoot ? 'bold' : 'normal')
        .attr('fill', '#1e293b')
        .attr('text-anchor', 'middle')
        .attr('dy', d => d.isRoot ? 28 : 20)
        .style('pointer-events', 'none')
        .text(d => d.name.length > 16 ? d.name.substring(0, 14) + '…' : d.name);

      simulation.on('tick', () => {
        link
          .attr('x1', d => d.source.x)
          .attr('y1', d => d.source.y)
          .attr('x2', d => d.target.x)
          .attr('y2', d => d.target.y);

        node.attr('cx', d => d.x).attr('cy', d => d.y);
        label.attr('x', d => d.x).attr('y', d => d.y);
      });

    } catch (err) {
      setError(err.message);
    } finally {
      setLoading(false);
    }
  }, [depth]);

  useEffect(() => {
    if (rootDrug) buildGraph(rootDrug);
  }, [rootDrug, depth, buildGraph]);

  const reset = () => {
    if (simulationRef.current) simulationRef.current.stop();
    d3.select(svgRef.current).selectAll('*').remove();
    setRootDrug(null);
    setQuery('');
    setError('');
  };

  return (
    <div className="p-6 flex flex-col h-full">
      <div className="flex items-center justify-between mb-4">
        <h1 className="text-2xl font-bold text-slate-900">Graph Explorer</h1>
        <button onClick={reset} className="px-3 py-1.5 text-sm border rounded hover:bg-slate-50">Reset</button>
      </div>

      {error && (
        <div className="mb-4 p-3 bg-red-50 border border-red-200 text-red-700 rounded text-sm">{error}</div>
      )}

      <div className="flex items-center gap-4 mb-4">
        {/* Drug search */}
        <div className="relative w-72">
          <input
            type="text"
            placeholder="Search root drug…"
            className="w-full border border-slate-300 rounded-lg px-4 py-2 text-sm focus:outline-none focus:ring-2 focus:ring-blue-500"
            value={query}
            onChange={e => { setQuery(e.target.value); setShowSugg(true); }}
            onFocus={() => setShowSugg(true)}
            onBlur={() => setTimeout(() => setShowSugg(false), 150)}
          />
          {showSugg && suggestions.length > 0 && (
            <div className="absolute top-full left-0 right-0 bg-white border border-slate-200 rounded-lg shadow-lg mt-1 z-10 max-h-48 overflow-y-auto">
              {suggestions.map(d => (
                <button
                  key={d.id}
                  onMouseDown={() => { setRootDrug(d); setQuery(d.name); setShowSugg(false); }}
                  className="w-full text-left px-4 py-2 hover:bg-slate-50 text-sm"
                >
                  {d.name} <span className="text-slate-400 text-xs">{d.id}</span>
                </button>
              ))}
            </div>
          )}
        </div>

        {/* Depth toggle */}
        <div className="flex gap-1">
          {[1, 2].map(d => (
            <button
              key={d}
              onClick={() => setDepth(d)}
              className={`px-3 py-1.5 text-sm rounded border transition-colors ${
                depth === d ? 'bg-[#0f172a] text-white border-[#0f172a]' : 'text-slate-600 hover:bg-slate-50'
              }`}
            >
              {d}-hop
            </button>
          ))}
        </div>

        {loading && (
          <div className="w-5 h-5 border-2 border-slate-300 border-t-blue-600 rounded-full animate-spin" />
        )}

        {/* Legend */}
        <div className="ml-auto flex items-center gap-4 text-xs text-slate-600">
          {['Major', 'Moderate'].map(s => (
            <span key={s} className="flex items-center gap-1">
              <span className="w-3 h-3 rounded-full inline-block" style={{ background: SEVERITY_COLOR[s] }} />
              {s}
            </span>
          ))}
          <span className="flex items-center gap-1">
            <span className="w-4 h-4 rounded-full inline-block bg-[#0f172a]" />
            Root
          </span>
        </div>
      </div>

      <div className="flex-1 bg-white border border-slate-200 rounded-xl overflow-hidden relative" style={{ minHeight: 500 }}>
        <svg ref={svgRef} width="100%" height="100%" style={{ minHeight: 500 }} />
        {!rootDrug && !loading && (
          <div className="absolute inset-0 flex items-center justify-center text-slate-400 pointer-events-none">
            <p>Search for a drug to explore its interaction network</p>
          </div>
        )}
      </div>

      {/* Tooltip */}
      {tooltip && (
        <div
          className="fixed z-50 bg-slate-900 text-white text-xs px-3 py-2 rounded-lg shadow-lg pointer-events-none max-w-xs"
          style={{ left: tooltip.x + 12, top: tooltip.y - 30 }}
        >
          {tooltip.content}
        </div>
      )}
    </div>
  );
}
