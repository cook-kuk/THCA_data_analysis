/* Plotly theme — medical data palette, dark-first */
window.THYROID_PLOTLY_TEMPLATE = {
  layout: {
    paper_bgcolor: 'rgba(0,0,0,0)',
    plot_bgcolor:  'rgba(0,0,0,0)',
    font: { family: 'Inter, system-ui, sans-serif', color: '#e2e8f0', size: 12 },
    colorway: ['#14b8a6','#2dd4bf','#5eead4','#38bdf8','#a5b4fc','#f59e0b','#fb7185','#e879f9'],
    margin: { l: 56, r: 24, t: 48, b: 56 },
    xaxis: { gridcolor: 'rgba(148,163,184,0.14)', zerolinecolor: 'rgba(148,163,184,0.12)', linecolor: 'rgba(148,163,184,0.18)' },
    yaxis: { gridcolor: 'rgba(148,163,184,0.14)', zerolinecolor: 'rgba(148,163,184,0.12)', linecolor: 'rgba(148,163,184,0.18)' },
    legend: { bgcolor: 'rgba(11,26,58,0.5)', bordercolor: 'rgba(45,212,191,0.2)', borderwidth: 1 },
    hoverlabel: { bgcolor: 'rgba(7,19,43,.94)', bordercolor: 'rgba(45,212,191,.4)', font: { color: '#e2e8f0' } }
  }
};
