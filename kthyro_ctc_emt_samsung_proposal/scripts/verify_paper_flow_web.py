#!/usr/bin/env python3
from __future__ import annotations

import asyncio
from pathlib import Path

from playwright.async_api import async_playwright


URL = "http://40.82.129.113:8012/papers_hub_2026_05_04/kthyro_ctc_emt_paper_flow_master_dossier.html"
OUT = Path("/home/seungho/personal/THCA_data_analysis/kthyro_ctc_emt_samsung_proposal/outputs/paper_flow_master_dossier/figures")


async def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)
        page = await browser.new_page(viewport={"width": 1440, "height": 1100}, device_scale_factor=1)
        await page.goto(URL, wait_until="networkidle")
        await page.screenshot(path=str(OUT / "web_paper_flow_hero.png"), full_page=False)
        await page.locator("#figures").scroll_into_view_if_needed()
        await page.screenshot(path=str(OUT / "web_paper_flow_figures.png"), full_page=False)

        mobile = await browser.new_page(viewport={"width": 390, "height": 1100}, is_mobile=True)
        await mobile.goto(URL, wait_until="networkidle")
        await mobile.screenshot(path=str(OUT / "web_paper_flow_mobile.png"), full_page=False)
        await browser.close()

    print(OUT / "web_paper_flow_hero.png")
    print(OUT / "web_paper_flow_figures.png")
    print(OUT / "web_paper_flow_mobile.png")


if __name__ == "__main__":
    asyncio.run(main())
