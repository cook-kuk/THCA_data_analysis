/* Section 2 — D3 sunburst of cancer → classifier → fold. */
(function (global) {
  const WIDTH = 760;
  const HEIGHT = 760;
  const RADIUS = Math.min(WIDTH, HEIGHT) / 2 - 10;

  function dialColor(v) {
    if (v == null || isNaN(v)) return '#52525a';
    // green (0) → amber (0.3) → red (0.5)
    if (v < 0.05) return '#6AB04C';
    if (v < 0.1)  return '#A1B84A';
    if (v < 0.2)  return '#D4B045';
    if (v < 0.3)  return '#F5A623';
    if (v < 0.4)  return '#EF7E24';
    return '#e74c3c';
  }

  let svg, g, root, partition, arc, current, app, tooltip;

  function renderSunburst(containerEl, tree) {
    svg = d3.select(containerEl).append('svg')
      .attr('viewBox', `0 0 ${WIDTH} ${HEIGHT}`)
      .attr('width', '100%')
      .style('max-height', '80vh')
      .style('display', 'block');

    g = svg.append('g').attr('transform', `translate(${WIDTH/2},${HEIGHT/2})`);

    tooltip = d3.select('body').append('div').attr('class', 'sunburst-tooltip');

    root = d3.hierarchy(tree)
      .sum(d => d.children ? 0 : 1)
      .sort((a, b) => (b.value || 0) - (a.value || 0));

    partition = d3.partition().size([2 * Math.PI, root.height + 1])(root);

    arc = d3.arc()
      .startAngle(d => d.x0)
      .endAngle(d => d.x1)
      .padAngle(d => Math.min((d.x1 - d.x0) / 2, 0.004))
      .padRadius(RADIUS * 1.5)
      .innerRadius(d => d.y0 * (RADIUS / root.height))
      .outerRadius(d => Math.max(d.y0 * (RADIUS/root.height), d.y1 * (RADIUS/root.height) - 1));

    root.each(d => d.current = d);

    const paths = g.selectAll('path')
      .data(root.descendants().slice(1))
      .join('path')
      .attr('fill', d => {
        if (d.depth === 1) {
          return dialColor(d.data.dial_max);
        }
        if (d.depth === 2) return dialColor(d.data.dial);
        if (d.depth === 3) return dialColor(d.data.dial);
        return '#333';
      })
      .attr('fill-opacity', d => d.depth === 1 ? 0.9 : (d.depth === 2 ? 0.75 : 0.5))
      .attr('d', arc)
      .style('cursor', 'pointer')
      .on('mouseenter', (e, d) => {
        tooltip.classed('visible', true).html(formatTooltip(d));
      })
      .on('mousemove', (e) => {
        tooltip.style('left', (e.pageX + 14) + 'px').style('top', (e.pageY + 8) + 'px');
      })
      .on('mouseleave', () => tooltip.classed('visible', false))
      .on('click', (e, d) => {
        updateSelection(d);
        if (d.children) clicked(d);
      });

    // center label
    const label = g.append('text')
      .attr('text-anchor', 'middle')
      .attr('dy', '0.35em')
      .style('font-family', 'var(--font-mono)')
      .style('fill', 'var(--accent)')
      .style('font-size', '0.95rem')
      .text('DIAL specificity');

    const sub = g.append('text')
      .attr('text-anchor', 'middle')
      .attr('y', 20)
      .style('font-family', 'var(--font-mono)')
      .style('fill', 'var(--fg-2)')
      .style('font-size', '0.72rem')
      .text('click segment · hover for details');

    current = root;
    updateBreadcrumb(root);

    function clicked(p) {
      current = p;
      updateBreadcrumb(p);
      root.each(d => d.target = {
        x0: Math.max(0, Math.min(1, (d.x0 - p.x0) / (p.x1 - p.x0))) * 2 * Math.PI,
        x1: Math.max(0, Math.min(1, (d.x1 - p.x0) / (p.x1 - p.x0))) * 2 * Math.PI,
        y0: Math.max(0, d.y0 - p.depth),
        y1: Math.max(0, d.y1 - p.depth),
      });

      const tx = g.transition().duration(550);
      paths.transition(tx)
        .tween('data', d => {
          const i = d3.interpolate(d.current, d.target);
          return (tt) => { d.current = i(tt); };
        })
        .attrTween('d', d => () => arc(d.current));
    }
  }

  function formatTooltip(d) {
    const lines = [];
    if (d.depth === 1) {
      lines.push(`<div><span class="key">cancer:</span> <span class="val">${d.data.name}</span></div>`);
      lines.push(`<div><span class="key">max DIAL:</span> ${fmt(d.data.dial_max, 3)}</div>`);
      lines.push(`<div><span class="key">n samples:</span> ${d.data.n_samples || '—'}</div>`);
    } else if (d.depth === 2) {
      lines.push(`<div><span class="key">classifier:</span> <span class="val">${d.data.name}</span></div>`);
      lines.push(`<div><span class="key">DIAL:</span> ${fmt(d.data.dial, 3)}</div>`);
      lines.push(`<div><span class="key">AUC<sub>pre</sub>:</span> ${fmt(d.data.auc_pre, 3)}</div>`);
      lines.push(`<div><span class="key">AUC<sub>post</sub>:</span> ${fmt(d.data.auc_post, 3)}</div>`);
      lines.push(`<div><span class="key">interp:</span> ${d.data.interpretation}</div>`);
      lines.push(`<div><span class="key">flip:</span> ${d.data.is_label_flip ? 'yes' : 'no'}</div>`);
    } else if (d.depth === 3) {
      lines.push(`<div><span class="key">fold:</span> <span class="val">${d.data.name}</span></div>`);
      lines.push(`<div><span class="key">DIAL:</span> ${fmt(d.data.dial, 3)}</div>`);
      lines.push(`<div><span class="key">AUC<sub>post</sub>:</span> ${fmt(d.data.auc_post, 3)}</div>`);
      lines.push(`<div class="mt-1 opacity-70 text-xs">per-fold values estimated from aggregate</div>`);
    }
    return lines.join('');
  }

  function fmt(v, d = 3) {
    if (v == null || isNaN(v)) return '—';
    return Number(v).toFixed(d);
  }

  function updateBreadcrumb(node) {
    const crumbEl = document.getElementById('sunburst-breadcrumb');
    if (!crumbEl) return;
    const ancestors = node.ancestors().reverse();
    crumbEl.innerHTML = ancestors.map((a, i) => {
      const isLast = i === ancestors.length - 1;
      return `<span class="crumb ${isLast ? 'current' : ''}" data-depth="${a.depth}">${a.data.name}</span>` +
             (isLast ? '' : '<span class="sep">›</span>');
    }).join('');
    crumbEl.querySelectorAll('.crumb').forEach((el, i) => {
      el.addEventListener('click', () => {
        if (i === 0) {
          d3.select('#sunburst-container svg g').selectAll('path').each(function(d) {});
        }
      });
    });
  }

  function updateSelection(d) {
    // write details to the right panel
    let cancer = '—', clf = '—', dial = '—', ap = '—', app_ = '—', interp = '—', flip = '—', n = '—';
    if (d.depth === 1) {
      cancer = d.data.name;
      dial = fmt(d.data.dial_max, 3);
      n = d.data.n_samples || '—';
    } else if (d.depth === 2) {
      cancer = d.parent.data.name;
      clf = d.data.name;
      dial = fmt(d.data.dial, 3);
      ap = fmt(d.data.auc_pre, 3);
      app_ = fmt(d.data.auc_post, 3);
      interp = d.data.interpretation;
      flip = d.data.is_label_flip ? 'yes' : 'no';
    } else if (d.depth === 3) {
      cancer = d.parent.parent.data.name;
      clf = d.parent.data.name;
      dial = fmt(d.data.dial, 3);
      app_ = fmt(d.data.auc_post, 3);
    }
    setTxt('det-cancer', cancer);
    setTxt('det-clf', clf);
    setTxt('det-dial', dial);
    setTxt('det-aucpre', ap);
    setTxt('det-aucpost', app_);
    setTxt('det-interp', interp);
    setTxt('det-flip', flip);
    setTxt('det-n', n);
    setTxt('det-view', cancer);

    if (window.__v9) {
      window.__v9.sunburstSel = { cancer, classifier: clf };
    }
  }

  function setTxt(id, v) {
    const el = document.getElementById(id);
    if (el) el.textContent = String(v);
  }

  const Section2 = {
    init(_app) {
      app = _app;
      const el = document.getElementById('sunburst-container');
      if (!el || !app.data) return;
      renderSunburst(el, app.data.hierarchy_tree);
    },
  };
  global.Section2 = Section2;
})(window);
