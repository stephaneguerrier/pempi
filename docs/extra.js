// Promote each row of the "Where to Next" tile-grid table to a fully
// clickable card. The pkgdown extra.css selector
//   main h2 + table:not(:has(thead tr th:not(:empty)))
// renders certain markdown tables as a tile grid. Their <tr>s contain a
// link inside <td:nth-child(2)>; without this shim, only the link text is
// clickable, which makes the tile feel half-broken when users click the
// icon area or the description.

(function () {
  function isTileGridTable(table) {
    if (!table.previousElementSibling || table.previousElementSibling.tagName !== "H2") {
      return false;
    }
    var thead = table.querySelector("thead");
    if (!thead) {
      return true;
    }
    var ths = thead.querySelectorAll("th");
    for (var i = 0; i < ths.length; i++) {
      if (ths[i].textContent.trim().length > 0) {
        return false;
      }
    }
    return true;
  }

  function makeRowClickable(tr) {
    var link = tr.querySelector("a[href]");
    if (!link) {
      return;
    }
    var href = link.getAttribute("href");
    var target = link.getAttribute("target");
    var label = link.textContent.trim();

    tr.setAttribute("role", "link");
    tr.setAttribute("tabindex", "0");
    if (label) {
      tr.setAttribute("aria-label", label);
    }
    tr.style.cursor = "pointer";

    function navigate(openInNewTab) {
      if (openInNewTab || target === "_blank") {
        window.open(href, "_blank", "noopener");
      } else {
        window.location.href = href;
      }
    }

    tr.addEventListener("click", function (e) {
      if (e.target.closest("a")) {
        return;
      }
      if (window.getSelection && window.getSelection().toString().length > 0) {
        return;
      }
      navigate(e.metaKey || e.ctrlKey);
    });

    tr.addEventListener("keydown", function (e) {
      if (e.key === "Enter" || e.key === " ") {
        e.preventDefault();
        navigate(false);
      }
    });
  }

  function init() {
    var tables = document.querySelectorAll("main h2 + table");
    tables.forEach(function (table) {
      if (!isTileGridTable(table)) {
        return;
      }
      table.querySelectorAll("tbody tr").forEach(makeRowClickable);
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();

// Pkgdown hard-codes the reference-index page title as "Function reference".
// Override to "Function Reference" so the casing matches the home-page H1
// and the navbar entry.
(function () {
  function fixReferenceIndexTitle() {
    if (document.title === "Function reference • pempi") {
      document.title = "Function Reference • pempi";
    }
    document.querySelectorAll("h1").forEach(function (h1) {
      if (h1.textContent.trim() === "Function reference") {
        h1.textContent = "Function Reference";
      }
    });
  }
  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", fixReferenceIndexTitle);
  } else {
    fixReferenceIndexTitle();
  }
})();
