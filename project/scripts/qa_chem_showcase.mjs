/* Playwright QA for pages/platform_chem_showcase.html
   - tests three viewports: 390x844, 768x1024, 1440x900
   - captures console errors, network failures, horizontal overflow
   - saves screenshots into logs/chem_showcase_screenshots/
*/
import { chromium } from "playwright";
import fs from "node:fs";
import path from "node:path";

const URL = "http://127.0.0.1:8012/reports/html/pages/platform_chem_showcase.html";
const OUT = "/opt/thyroid-dash/project/logs/chem_showcase_screenshots";
fs.mkdirSync(OUT, { recursive: true });

const viewports = [
  { name: "mobile_390x844",   width: 390,  height: 844,  device: true },
  { name: "tablet_768x1024",  width: 768,  height: 1024, device: false },
  { name: "desktop_1440x900", width: 1440, height: 900,  device: false },
];

const summary = [];

const browser = await chromium.launch({ args: ["--no-sandbox"] });
try {
  for (const vp of viewports) {
    const ctx = await browser.newContext({
      viewport: { width: vp.width, height: vp.height },
      deviceScaleFactor: 1,
    });
    const page = await ctx.newPage();

    const consoleErrors = [];
    const failedRequests = [];
    page.on("console", msg => { if (msg.type() === "error") consoleErrors.push(msg.text()); });
    page.on("pageerror", err => consoleErrors.push("pageerror: " + err.message));
    page.on("requestfailed", req => failedRequests.push(req.url() + " -> " + req.failure()?.errorText));

    const t0 = Date.now();
    await page.goto(URL, { waitUntil: "domcontentloaded", timeout: 20000 });
    // wait for the main UI state (target dropdown populated)
    await page.waitForFunction(() => {
      const sel = document.querySelector("#target-select");
      return sel && sel.options && sel.options.length >= 8;
    }, { timeout: 10000 }).catch(() => {});
    // give charts/viewer time to mount
    await page.waitForTimeout(3000);
    const loadMs = Date.now() - t0;

    // Check horizontal overflow
    const horizOverflow = await page.evaluate(() => {
      const docWidth = document.documentElement.scrollWidth;
      const viewWidth = window.innerWidth;
      return { docWidth, viewWidth, overflow: docWidth - viewWidth };
    });

    // Check viewer initialized
    const viewerReady = await page.evaluate(() => {
      return !!(window.MolViewer && window.$3Dmol);
    });

    // Test target dropdown change
    let dropdownOk = false;
    try {
      await page.selectOption("#target-select", "LDLR");
      await page.waitForTimeout(1500);
      dropdownOk = await page.evaluate(() => window.viewerState && window.viewerState.currentTarget === "LDLR");
    } catch (e) { dropdownOk = false; }

    // Test style toggle
    let styleOk = false;
    try {
      await page.click('.style-toggle button[data-style="stick"]');
      await page.waitForTimeout(500);
      styleOk = await page.evaluate(() => window.viewerState.style === "stick");
    } catch (e) { styleOk = false; }

    // restore CYP1B1 for screenshot consistency
    try {
      await page.selectOption("#target-select", "CYP1B1");
      await page.click('.style-toggle button[data-style="cartoon"]');
      await page.waitForTimeout(1500);
    } catch (e) {}

    const shotPath = path.join(OUT, `${vp.name}.png`);
    await page.screenshot({ path: shotPath, fullPage: true });

    // Screenshot the viewer section specifically too
    const viewerBox = await page.locator("#binding-viewer").boundingBox();
    if (viewerBox) {
      await page.screenshot({
        path: path.join(OUT, `${vp.name}_viewer.png`),
        clip: {
          x: Math.max(0, viewerBox.x),
          y: Math.max(0, viewerBox.y),
          width: Math.min(vp.width, viewerBox.width),
          height: Math.min(1200, viewerBox.height),
        },
      });
    }

    summary.push({
      viewport: vp.name,
      load_ms: loadMs,
      viewer_ready: viewerReady,
      dropdown_ok: dropdownOk,
      style_toggle_ok: styleOk,
      horizontal_overflow_px: horizOverflow.overflow,
      doc_width: horizOverflow.docWidth,
      view_width: horizOverflow.viewWidth,
      console_errors: consoleErrors,
      failed_requests: failedRequests,
      screenshot: shotPath,
    });
    await ctx.close();
  }
} finally {
  await browser.close();
}

fs.writeFileSync(path.join(OUT, "qa_summary.json"), JSON.stringify(summary, null, 2));
console.log(JSON.stringify(summary, null, 2));
