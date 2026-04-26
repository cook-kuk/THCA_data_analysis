/* Section 5 — floating AI panel. If user provides Anthropic API key, stream
   real answers. Otherwise fall back to static FAQ answers. */
(function (global) {
  let app = null;
  const KEY_STORAGE = 'v9_anthropic_key';

  const FAQ = [
    {
      match: /(v5\.?1|leak|leakage|artifact|data leak|why.*0\.494)/i,
      answer: `<b>v5.1 leakage:</b> ComBat was fit on the full pooled X (train + test cohorts) before
               the LODO split. Even though the classifier never saw the held-out cohort during training,
               the input matrix had already been adjusted with statistics that included the held-out
               samples. The empirical-Bayes priors (α, β, γ, δ) borrowed strength across all rows. So
               LODO appeared to invert the THCA classifier (AUC_post = 0.006, DIAL = 0.494) — but this was
               leakage, not biology. Under v5.2 (per-fold ComBat) THCA AUC_post = 0.995, DIAL = 0.000.`,
    },
    {
      match: /(why|reason).*(thca|thyroid).*(flip|invert)/i,
      answer: `<b>It doesn't (under v5.2).</b> v5.1 reported THCA DIAL = 0.494 in 4/5 classifiers, but
               that "flip" was a ComBat-on-pooled-X artifact. With proper per-fold ComBat (fit on training
               cohort only, held-out cohort centered locally without seeing labels), THCA gives DIAL = 0.000
               and AUC_post = 0.995 in all 5 classifiers. ΔDIAL = −0.494 across the board.`,
    },
    {
      match: /(why|reason).*(lgg|glioma|skcm|luad|coad).*(not flip|doesn't flip|true biology|unchanged)/i,
      answer: `LGG (IDHmut vs IDHwt) has strong, monogenic transcriptional consequences. SKCM/LUAD/COAD
               also returned no flip in <em>both</em> v5.1 and v5.2 — they were never affected by the
               leak because their cohort statistics happened to align with their label structure. THCA
               was the only cancer where the v5.1 leak produced an apparent flip, and it was exactly
               there that v5.2 corrected the call.`,
    },
    {
      match: /(should|still|use|valid).*dial/i,
      answer: `<b>Yes — but its role has changed.</b> Under the v5.1 protocol DIAL was a measurement
               instrument; under the v5.2 reframing it's a leakage detector. If your pipeline reports
               high DIAL, the most parsimonious explanation is that ComBat (or another corrector) is
               seeing the test cohort. The v5.2 paper repositions DIAL as a methodological audit tool,
               not a biology metric.`,
    },
    {
      match: /(combat-seq|pycombat|reference).*(difference|instead|vs)/i,
      answer: `ComBat-seq expects integer counts (NB distribution); ComBat (pycombat) assumes Gaussian
               log-expression. Both have the same leakage pitfall: if you fit on the full matrix before
               splitting, you leak. Use per-fold fitting regardless of estimator family.`,
    },
    {
      match: /(which|what).*gene.*(drive|flip|thca)/i,
      answer: `Top 3 TCGA-vs-GEO genes in THCA: TG, CTSB, FN1 — TCGA−GEO mean ≥ 8 log2 units. Under v5.1
               these features drove the apparent flip. Under v5.2 they're correctly removed from the
               training matrix only; the held-out cohort is locally centered. The signal that matters
               (BRAF-vs-RAS) is preserved.`,
    },
    {
      match: /(interpret|what is).*dial/i,
      answer: `DIAL = Detected Inversion of Axis after Label-preserving correction. Formula:
               max(0, (1 − AUC_post) − 0.5) when AUC_post < 0.5. v5.1 reported DIAL > 0.3 ("batch_entangled")
               on THCA. v5.2 self-audit: that DIAL was a leakage signal — the metric correctly fired on a
               leaky pipeline, which is a methodological win, not a biological finding.`,
    },
  ];

  function appendMessage(role, html) {
    const body = document.getElementById('ai-body');
    const d = document.createElement('div');
    d.className = 'ai-message ' + role;
    d.innerHTML = html;
    body.appendChild(d);
    body.scrollTop = body.scrollHeight;
  }

  function staticAnswer(q) {
    for (const f of FAQ) {
      if (f.match.test(q)) return f.answer;
    }
    return `I don't have a pre-written answer for that. Paste your Anthropic API key (🔑 button) for a
            live answer. The data I have loaded: ${(app.data?.dial_results || []).length} DIAL rows across
            ${(app.data?.meta?.cancers || []).join(', ')}, and 20 batch-discriminant genes per cancer.`;
  }

  async function callClaude(question, apiKey) {
    const summary = makeContextSummary();
    const body = {
      model: 'claude-sonnet-4-6',
      max_tokens: 600,
      messages: [{
        role: 'user',
        content:
`Context — v5.1 DIAL results:
${summary}

Question: ${question}

Answer concisely in plain prose. Cite specific numbers from the context if relevant.`,
      }],
    };

    try {
      const res = await fetch('https://api.anthropic.com/v1/messages', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'x-api-key': apiKey,
          'anthropic-version': '2023-06-01',
          'anthropic-dangerous-direct-browser-access': 'true',
        },
        body: JSON.stringify(body),
      });
      if (!res.ok) {
        const err = await res.text();
        return `<span class="text-red-400">API error (${res.status}): ${err.slice(0, 300)}</span>`;
      }
      const j = await res.json();
      return (j.content?.[0]?.text || '(empty response)').replace(/\n/g, '<br>');
    } catch (e) {
      return `<span class="text-red-400">Network error: ${e.message}</span>`;
    }
  }

  function makeContextSummary() {
    if (!app.data) return '(data not loaded)';
    const rows = app.data.dial_results;
    const lines = rows.slice(0, 25).map(r =>
      `- ${r.cancer}/${r.classifier}: DIAL=${(r.dial||0).toFixed(3)}, AUC_pre=${(r.auc_pre||0).toFixed(3)}, AUC_post=${(r.auc_post||0).toFixed(3)}, ${r.interpretation}`);
    return lines.join('\n');
  }

  function promptForKey() {
    const existing = localStorage.getItem(KEY_STORAGE) || '';
    const k = prompt('Paste your Anthropic API key (stored in this browser only). Leave empty to clear.', existing);
    if (k === null) return;
    if (k === '') {
      localStorage.removeItem(KEY_STORAGE);
      appendMessage('bot', 'API key cleared. Falling back to static FAQ answers.');
    } else {
      localStorage.setItem(KEY_STORAGE, k);
      appendMessage('bot', 'API key saved in browser. Ask away.');
    }
  }

  async function handleSend() {
    const input = document.getElementById('ai-input');
    const q = input.value.trim();
    if (!q) return;
    input.value = '';
    appendMessage('user', escapeHtml(q));
    const key = localStorage.getItem(KEY_STORAGE);
    if (key) {
      appendMessage('bot', '<span class="opacity-60">thinking…</span>');
      const ans = await callClaude(q, key);
      const last = document.querySelector('.ai-message.bot:last-child');
      if (last) last.innerHTML = ans;
    } else {
      appendMessage('bot', staticAnswer(q));
    }
  }

  function escapeHtml(s) {
    return s.replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' })[c]);
  }

  const Section5 = {
    init(_app) {
      app = _app;
      document.getElementById('ai-send')?.addEventListener('click', handleSend);
      document.getElementById('ai-input')?.addEventListener('keydown', e => {
        if (e.key === 'Enter' && !e.shiftKey) { e.preventDefault(); handleSend(); }
      });
      document.getElementById('ai-key-btn')?.addEventListener('click', promptForKey);
    },
  };
  global.Section5 = Section5;
})(window);
