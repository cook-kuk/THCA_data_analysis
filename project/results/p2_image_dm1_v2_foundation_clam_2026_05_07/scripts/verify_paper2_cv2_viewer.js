const { chromium } = require("playwright");
const path = require("path");

(async () => {
  const root = "/home/seungho/personal/THCA_data_analysis";
  const html = path.join(root, "project/papers_hub_2026_05_04/paper2_cv2_visual_summary.html");
  const out = path.join(
    root,
    "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/analysis_supp/paper2_cv2_visual_summary_2026_05_10/paper2_cv2_interactive_viewer_verify.png"
  );
  const browser = await chromium.launch({ headless: true });
  const page = await browser.newPage({ viewport: { width: 1440, height: 1050 }, deviceScaleFactor: 1 });
  await page.goto(`file://${html}`, { waitUntil: "networkidle" });
  await page.locator("#figure-explorer").scrollIntoViewIfNeeded();
  await page.locator(".figure-option").nth(3).click();
  await page.locator("[data-zoom-in]").click();
  await page.locator("[data-zoom-in]").click();
  const title = await page.locator("#figureTitle").textContent();
  const width = await page.locator("#figureImage").evaluate((img) => getComputedStyle(img).width);
  await page.locator("#figure-explorer").screenshot({ path: out });
  await browser.close();
  console.log(JSON.stringify({ ok: true, title, width, screenshot: out }, null, 2));
})();
