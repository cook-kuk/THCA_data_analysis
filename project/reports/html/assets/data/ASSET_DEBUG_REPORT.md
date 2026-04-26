# ASSET_DEBUG_REPORT

## broken paths found
- index homepage nav first link used ../index.html in generated source
- pages command index used pages/... paths even inside pages/ context

## CDN replacements
- tailwind -> local vendor script
- alpine -> local vendor script
- aos css/js -> local vendor copies
- plotly -> local vendor copy
- marked -> local vendor copy

## local vendor fallbacks created
- project/reports/html/assets/vendor/plotly/
- project/reports/html/assets/vendor/tailwind/
- project/reports/html/assets/vendor/alpine/
- project/reports/html/assets/vendor/aos/
- project/reports/html/assets/vendor/marked/
- project/reports/html/assets/vendor/fuse/
- project/reports/html/assets/vendor/lucide/

## cache/CORS/path fixes
- viewport/meta refreshed in output html
- HTML now points to local vendor assets for key libs
- runtime path normalization added in main.js
