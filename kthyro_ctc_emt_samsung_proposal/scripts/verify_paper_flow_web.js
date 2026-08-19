const { chromium } = require('playwright');
const path = require('path');

const url = 'http://40.82.129.113:8012/papers_hub_2026_05_04/kthyro_ctc_emt_paper_flow_master_dossier.html';
const out = '/home/seungho/personal/THCA_data_analysis/kthyro_ctc_emt_samsung_proposal/outputs/paper_flow_master_dossier/figures';

(async () => {
  const browser = await chromium.launch({ headless: true });
  const page = await browser.newPage({ viewport: { width: 1440, height: 1100 }, deviceScaleFactor: 1 });
  await page.goto(url, { waitUntil: 'networkidle' });
  await page.screenshot({ path: path.join(out, 'web_paper_flow_hero.png'), fullPage: false });
  await page.locator('#figures').scrollIntoViewIfNeeded();
  await page.screenshot({ path: path.join(out, 'web_paper_flow_figures.png'), fullPage: false });
  const mobile = await browser.newPage({ viewport: { width: 390, height: 1100 }, isMobile: true });
  await mobile.goto(url, { waitUntil: 'networkidle' });
  await mobile.screenshot({ path: path.join(out, 'web_paper_flow_mobile.png'), fullPage: false });
  await browser.close();
  console.log(path.join(out, 'web_paper_flow_hero.png'));
  console.log(path.join(out, 'web_paper_flow_figures.png'));
  console.log(path.join(out, 'web_paper_flow_mobile.png'));
})();
