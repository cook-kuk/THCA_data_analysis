window.ThyroidDash = window.ThyroidDash || {};
window.ThyroidDash.filterPanel = {
  save(state) {
    localStorage.setItem("thyroid-dash-filters", JSON.stringify(state || {}));
  },
  load() {
    try { return JSON.parse(localStorage.getItem("thyroid-dash-filters") || "{}"); }
    catch (e) { return {}; }
  },
  share() {
    navigator.clipboard.writeText(location.href);
  }
};
