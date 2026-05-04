import React from 'react';
import { BrowserRouter, Routes, Route, NavLink } from 'react-router-dom';
import DrugSearch from './pages/DrugSearch.jsx';
import InteractionChecker from './pages/InteractionChecker.jsx';
import GraphExplorer from './pages/GraphExplorer.jsx';
import CRUDManager from './pages/CRUDManager.jsx';
import Analytics from './pages/Analytics.jsx';

const navItems = [
  { to: '/', label: 'Drug Search', end: true },
  { to: '/checker', label: 'Interaction Checker' },
  { to: '/graph', label: 'Graph Explorer' },
  { to: '/manage', label: 'CRUD Manager' },
  { to: '/analytics', label: 'Analytics' },
];

export default function App() {
  return (
    <BrowserRouter>
      <div className="flex min-h-screen bg-slate-50">
        {/* Sidebar */}
        <aside className="w-56 bg-[#0f172a] text-white flex flex-col shrink-0">
          <div className="px-5 py-6 border-b border-slate-700">
            <h1 className="text-lg font-bold tracking-tight">DDI Explorer</h1>
            <p className="text-xs text-slate-400 mt-1">Drug Interaction Graph</p>
          </div>
          <nav className="flex-1 py-4">
            {navItems.map(({ to, label, end }) => (
              <NavLink
                key={to}
                to={to}
                end={end}
                className={({ isActive }) =>
                  `block px-5 py-2.5 text-sm transition-colors ${
                    isActive
                      ? 'bg-slate-700 text-white font-medium'
                      : 'text-slate-400 hover:text-white hover:bg-slate-800'
                  }`
                }
              >
                {label}
              </NavLink>
            ))}
          </nav>
        </aside>

        {/* Main content */}
        <main className="flex-1 overflow-auto">
          <Routes>
            <Route path="/" element={<DrugSearch />} />
            <Route path="/checker" element={<InteractionChecker />} />
            <Route path="/graph" element={<GraphExplorer />} />
            <Route path="/manage" element={<CRUDManager />} />
            <Route path="/analytics" element={<Analytics />} />
          </Routes>
        </main>
      </div>
    </BrowserRouter>
  );
}
