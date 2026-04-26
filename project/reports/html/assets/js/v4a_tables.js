/* v4a_tables.js - render TSV tables with numeric formatting + show-all toggle.
   Not affiliated with any drug-discovery AI company.
*/
(function(){
  'use strict';

  function isNumeric(v){
    if (v === '' || v == null) return false;
    // allow scientific, signed, decimal
    return /^-?\d+(\.\d+)?([eE][-+]?\d+)?$/.test(String(v).trim());
  }

  function formatNumber(v){
    if (v === '' || v == null) return '';
    var n = Number(v);
    if (!isFinite(n)) return v;
    var a = Math.abs(n);
    if (n === 0) return '0';
    if (a > 0 && a < 1e-2) return n.toExponential(2);
    if (a >= 1e4) return n.toExponential(2);
    if (Number.isInteger(n) && a < 1e6) return String(n);
    // 3 significant figures
    return n.toPrecision(3).replace(/\.?0+$/, '');
  }

  function parseTSV(text){
    var lines = text.replace(/\r/g, '').split('\n').filter(function(l){ return l.length > 0; });
    if (!lines.length) return { header: [], rows: [] };
    var header = lines[0].split('\t');
    var rows = [];
    for (var i = 1; i < lines.length; i++){
      rows.push(lines[i].split('\t'));
    }
    return { header: header, rows: rows };
  }

  function detectNumericCols(rows, header){
    var n = header.length;
    var numeric = new Array(n).fill(true);
    for (var r = 0; r < Math.min(rows.length, 80); r++){
      var row = rows[r];
      for (var c = 0; c < n; c++){
        var v = row[c];
        if (v === undefined || v === '' || v == null) continue;
        if (!isNumeric(v)) numeric[c] = false;
      }
    }
    return numeric;
  }

  function renderTable(container, header, rows, numericCols, topN, label){
    var hasMore = rows.length > topN;
    var shown = hasMore ? rows.slice(0, topN) : rows;

    var wrap = document.createElement('div');
    wrap.className = 'v4a-table-wrap';
    var table = document.createElement('table');
    table.className = 'v4a-table';
    var thead = document.createElement('thead');
    var trh = document.createElement('tr');
    header.forEach(function(h, i){
      var th = document.createElement('th');
      th.textContent = h;
      if (numericCols[i]) th.className = 'num';
      trh.appendChild(th);
    });
    thead.appendChild(trh);
    table.appendChild(thead);

    var tbody = document.createElement('tbody');
    function appendRows(dataRows){
      dataRows.forEach(function(r){
        var tr = document.createElement('tr');
        for (var i = 0; i < header.length; i++){
          var td = document.createElement('td');
          var v = r[i] === undefined ? '' : r[i];
          if (numericCols[i]){
            td.className = 'num';
            td.textContent = formatNumber(v);
          } else {
            td.textContent = v;
          }
          tr.appendChild(td);
        }
        tbody.appendChild(tr);
      });
    }
    appendRows(shown);
    table.appendChild(tbody);
    wrap.appendChild(table);
    container.appendChild(wrap);

    var footer = document.createElement('div');
    footer.className = 'v4a-table-footer';
    var span = document.createElement('span');
    span.textContent = (label || 'Rows') + ' · ' + shown.length + ' / ' + rows.length;
    footer.appendChild(span);
    if (hasMore){
      var btn = document.createElement('button');
      btn.type = 'button';
      btn.textContent = 'Show all (' + rows.length + ')';
      var expanded = false;
      btn.addEventListener('click', function(){
        if (!expanded){
          appendRows(rows.slice(topN));
          btn.textContent = 'Show top ' + topN;
          span.textContent = (label || 'Rows') + ' · ' + rows.length + ' / ' + rows.length;
          expanded = true;
        } else {
          // rebuild only top N
          tbody.innerHTML = '';
          appendRows(rows.slice(0, topN));
          btn.textContent = 'Show all (' + rows.length + ')';
          span.textContent = (label || 'Rows') + ' · ' + topN + ' / ' + rows.length;
          expanded = false;
        }
      });
      footer.appendChild(btn);
    }
    container.appendChild(footer);
  }

  function renderError(container, msg){
    var d = document.createElement('div');
    d.className = 'v4a-source-missing';
    d.innerHTML = '<strong>SOURCE UNAVAILABLE</strong>' + msg;
    container.appendChild(d);
  }

  window.renderTSVTable = function(containerId, tsvPath, opts){
    opts = opts || {};
    var topN = opts.topN || 30;
    var label = opts.label || 'rows';
    var sortBy = opts.sortBy; // column name to sort by desc (abs value)
    var container = document.getElementById(containerId);
    if (!container) return;
    fetch(tsvPath)
      .then(function(r){
        if (!r.ok) throw new Error('HTTP ' + r.status);
        return r.text();
      })
      .then(function(text){
        var parsed = parseTSV(text);
        if (!parsed.header.length){
          renderError(container, tsvPath + ' (empty file)');
          return;
        }
        var numeric = detectNumericCols(parsed.rows, parsed.header);
        var rows = parsed.rows;
        if (sortBy){
          var idx = parsed.header.indexOf(sortBy);
          if (idx >= 0){
            rows = rows.slice().sort(function(a, b){
              var av = Math.abs(parseFloat(a[idx])); if (!isFinite(av)) av = -Infinity;
              var bv = Math.abs(parseFloat(b[idx])); if (!isFinite(bv)) bv = -Infinity;
              return bv - av;
            });
          }
        }
        renderTable(container, parsed.header, rows, numeric, topN, label);
      })
      .catch(function(err){
        renderError(container, tsvPath + ' (' + err.message + ')');
      });
  };
})();
