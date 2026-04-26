const { chromium } = require('playwright');
const fs = require('fs');
const path = require('path');
(async()=>{
  const outDir = path.resolve('project/logs/ui_screenshots');
  fs.mkdirSync(outDir, { recursive: true });
  const browser = await chromium.launch({ headless: true });
  const results = [];
  for (const vp of [
    { name: 'desktop', width: 1440, height: 900 },
    { name: 'tablet', width: 768, height: 1024 },
    { name: 'mobile', width: 390, height: 844, isMobile: true }
  ]) {
    const context = await browser.newContext({ viewport: { width: vp.width, height: vp.height }, isMobile: !!vp.isMobile });
    const page = await context.newPage();
    const consoleErrors = [];
    const failed = [];
    page.on('console', msg => { if (msg.type() === 'error') consoleErrors.push(msg.text()); });
    page.on('requestfailed', req => failed.push(req.url()));
    await page.goto('http://127.0.0.1:8012/reports/html/index.html', { waitUntil: 'domcontentloaded', timeout: 60000 });
    await page.waitForTimeout(1200);
    const menu = page.locator('#td-menu-toggle');
    const menuCount = await menu.count();
    let drawerVisible = false;
    if (menuCount) {
      await menu.click();
      await page.waitForTimeout(250);
      drawerVisible = await page.locator('#td-mobile-drawer').isVisible().catch(()=>false);
      const backdrop = page.locator('#td-mobile-backdrop');
      if (await backdrop.count()) await backdrop.click({ position: { x: 10, y: 10 } }).catch(()=>{});
    }
    await page.screenshot({ path: path.join(outDir, `${vp.name}.png`), fullPage: true });
    const overflow = await page.evaluate(() => document.documentElement.scrollWidth > window.innerWidth + 2);
    const uiHealth = await page.evaluate(() => window.__uiHealth || null);
    results.push({ viewport: vp.name, menuCount, drawerVisible, overflow, uiHealth, consoleErrors, failed });
    await context.close();
  }
  await browser.close();
  fs.writeFileSync(path.resolve('project/logs/ui_playwright_results.json'), JSON.stringify(results, null, 2));
})();
