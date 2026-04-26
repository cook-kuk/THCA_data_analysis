// Thorough UI bug hunt across every dashboard page: console errors, runtime errors,
// failed network requests, render-time overflow. Skips npm deps — uses playwright already
// installed at the repo root node_modules/.
const { chromium } = require('playwright');
const fs = require('fs');
const path = require('path');

const BASE = process.env.BASE_URL || 'http://127.0.0.1:8012/reports/html';
const PAGES = [
  '/index.html',
  '/dashboard_v2_interactive.html',
  '/pages/01_overview.html',
  '/pages/02_datasets.html',
  '/pages/03_sample_master.html',
  '/pages/04_gene_panels.html',
  '/pages/05_eda.html',
  '/pages/06_scores.html',
  '/pages/07_ml_baseline.html',
  '/pages/08_panel_comparison.html',
  '/pages/09_shap.html',
  '/pages/10_gene_explorer.html',
  '/pages/11_cohort_compare.html',
  '/pages/12_business.html',
  '/pages/13_reports.html',
  '/pages/14_caveats.html',
];

(async () => {
  const browser = await chromium.launch({ headless: true, args: ['--no-sandbox'] });
  const context = await browser.newContext({ viewport: { width: 1440, height: 900 } });
  const results = [];
  for (const p of PAGES) {
    const page = await context.newPage();
    const consoleErrors = [];
    const pageErrors = [];
    const failed = [];
    page.on('console', msg => { if (msg.type() === 'error') consoleErrors.push(msg.text()); });
    page.on('pageerror', err => pageErrors.push(err.message || String(err)));
    page.on('requestfailed', req => failed.push({ url: req.url(), err: req.failure()?.errorText }));
    let status = 0;
    try {
      const resp = await page.goto(BASE + p, { waitUntil: 'domcontentloaded', timeout: 30000 });
      status = resp ? resp.status() : -1;
      await page.waitForTimeout(1500);
    } catch (e) {
      pageErrors.push('navigation: ' + e.message);
    }
    const overflow = await page.evaluate(() => document.documentElement.scrollWidth > window.innerWidth + 2).catch(()=>null);
    const h1Count = await page.evaluate(() => document.querySelectorAll('h1').length).catch(()=>null);
    const imgNoAlt = await page.evaluate(() => Array.from(document.querySelectorAll('img')).filter(i => !i.alt).length).catch(()=>null);
    // Count iframes with no src, EXCLUDING the figure-modal iframe whose src is intentionally
    // empty until the user clicks a figure card. (id="figure-modal-frame" is a known good pattern.)
    const brokenIframes = await page.evaluate(() => Array.from(document.querySelectorAll('iframe')).filter(f => !f.src && f.id !== 'figure-modal-frame').length).catch(()=>null);
    // Filter out benign console noise
    const localFailed = failed.filter(f => f.url && !f.url.startsWith('http://127.0.0.1') && !f.url.startsWith('http://localhost') ? false : true);
    results.push({ page: p, status, overflow, h1Count, imgNoAlt, brokenIframes, consoleErrors: consoleErrors.slice(0,6), pageErrors: pageErrors.slice(0,6), failed: failed.slice(0,6) });
    await page.close();
  }
  await context.close();
  await browser.close();
  const outPath = path.resolve(path.dirname(__filename), '..', 'logs', 'ui_bug_hunt.json');
  fs.mkdirSync(path.dirname(outPath), { recursive: true });
  fs.writeFileSync(outPath, JSON.stringify(results, null, 2));
  // Compact stdout report
  let totalErrors = 0, totalPageErrors = 0, totalFailed = 0;
  for (const r of results) { totalErrors += r.consoleErrors.length; totalPageErrors += r.pageErrors.length; totalFailed += r.failed.length; }
  console.log(`Scanned ${results.length} pages.  console-errors=${totalErrors}  page-errors=${totalPageErrors}  failed-requests=${totalFailed}`);
  for (const r of results) {
    const issues = [];
    if (r.status !== 200) issues.push(`status=${r.status}`);
    if (r.consoleErrors.length) issues.push(`console=${r.consoleErrors.length}`);
    if (r.pageErrors.length) issues.push(`pageerr=${r.pageErrors.length}`);
    if (r.failed.length) issues.push(`failed=${r.failed.length}`);
    if (r.overflow) issues.push('h-overflow');
    if (r.imgNoAlt) issues.push(`img-no-alt=${r.imgNoAlt}`);
    if (r.brokenIframes) issues.push(`empty-iframe=${r.brokenIframes}`);
    if (issues.length) {
      console.log(`  !! ${r.page}  ${issues.join(' ')}`);
      for (const e of r.consoleErrors.slice(0,3)) console.log('       console:', e.slice(0,180));
      for (const e of r.pageErrors.slice(0,3))   console.log('       pageerr:', e.slice(0,180));
      for (const f of r.failed.slice(0,3))       console.log('       failed :', f.url, '->', f.err);
    } else {
      console.log(`  OK ${r.page}`);
    }
  }
  console.log('Detail saved to:', outPath);
})();
