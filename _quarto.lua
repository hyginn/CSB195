-- _quarto.lua

-- Place functions that should execute every time a quarto document
-- is rendered in this project into this file. If you want to include a whole
-- directory of functions, you can require() it by placing code here, like:
-- local util = require('./src/lua/utils')

-- copyBox function
-- Transform fenced code blocks {.copyBox label="…"} into a styled container
-- with a copy button that copies only the <pre> payload.
-- Usage in .qmd:
-- ```{.copyBox label="<your label>t"}
-- Your literal text...
-- ```

-- -------- Helpers (top-level, local) --------

local function is_html()
  -- Prefer Quarto's helper when available; fallback to Pandoc's FORMAT
  if quarto and quarto.doc then
    return quarto.doc.is_format("html")
  end
  return FORMAT and FORMAT:match("html")
end

local function html_escape(s)
  -- Escape raw characters for safe inclusion inside <pre>…</pre>
  return (s:gsub("&", "&amp;")
           :gsub("<", "&lt;")
           :gsub(">", "&gt;"))
end

local function has_class(el, cls)
  local classes = el.classes or (el.attr and el.attr.classes) or {}
  -- Pandoc's List often has :includes; otherwise, fall back to a manual scan
  if classes.includes and classes:includes(cls) then return true end
  for _, c in ipairs(classes) do
    if c == cls then return true end
  end
  return false
end

local function get_attr(el, key)
  local attrs = el.attributes or (el.attr and el.attr.attributes) or {}
  return attrs[key]
end

-- -------- Core transformer --------

local function handle_copybox(el)
  if not is_html() then return nil end
  if not has_class(el, "copyBox") then return nil end

  local label   = get_attr(el, "label") or "Copy"
  local payload = el.text or ""  -- CodeBlock.text is a STRING

  local button = [[
<button class="copyButton" type="button" aria-label="Copy" onclick="copyContent(this)">
  <i class="bi bi-clipboard"></i>
</button>
]]

  local html = string.format([[
<div class="copyBox" label="%s">
  %s
  <pre>%s</pre>
</div>]],
    pandoc.utils.stringify(label),
    button,
    html_escape(payload)
  )

  return pandoc.RawBlock("html", html)
end

-- -------- Expose handlers --------

return {
  CodeBlock = function(el)
    return handle_copybox(el)  -- return nil for non-matching blocks
  end
}

-- [END]
