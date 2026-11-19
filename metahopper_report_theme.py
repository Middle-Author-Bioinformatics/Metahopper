"""Report palette for standalone and hosted light/dark views."""

THEME_CSS = '''
:root{--surface:#f4f7f8;--hover:#f0f6f7;--grid:#e3eaed;--on-acc:#fff;--ok-bg:#e1f1e9;--warn-bg:#f9e5e2;--review:#8b5d16;--review-bg:#fff0d9;--pseudo:#a32626;--pseudo-bg:#fce8e6;--selection:#e2f0f2}
:root[data-theme="dark"]{color-scheme:dark;--bg:#10191e;--card:#18252d;--surface:#1e2e37;--ink:#e2edf1;--mut:#a2b6c1;--line:#344953;--acc:#79d1d7;--accs:#253d47;--hover:#20343e;--grid:#344953;--on-acc:#102c32;--ok:#9bdbb7;--ok-bg:#203f34;--warn:#ffaca3;--warn-bg:#482c2d;--review:#f2cb83;--review-bg:#443a25;--pseudo:#ffaaa6;--pseudo-bg:#492e32;--selection:#284651}
:root a{color:var(--acc)}
:root body{background:var(--bg);color:var(--ink)}
:root tbody tr:hover{background:var(--hover)}
:root dialog,:root canvas,:root .graph-stage{background:var(--card);color:var(--ink)}
:root button.p,:root .tab.on{background:var(--acc);color:var(--on-acc);border-color:var(--acc)}
:root button.d,:root .asset-link{background:var(--card)}
:root .pill.ok,:root .feature-tag.named,:root .pathway-tag.ok{background:var(--ok-bg);color:var(--ok)}
:root .pill.no,:root .pathway-tag.bad{background:var(--warn-bg);color:var(--warn)}
:root .feature-tag,:root .pathway-tag{background:var(--accs);color:var(--mut)}
:root .feature-tag.hypothetical,:root .pathway-tag.partial,:root .pathway-tag.missing{background:var(--review-bg);color:var(--review)}
:root .feature-tag.pseudogene{background:var(--pseudo-bg);color:var(--pseudo)}
:root .feature-pseudogene .annotation-feature-link,:root .feature-pseudogene .feature-gene{color:var(--pseudo)}
:root .sequence-preview,:root .pathway-sequence{background:var(--surface);color:var(--ink)}
:root .pathway-warning{background:var(--review-bg);color:var(--ink)}
:root .figure-host line{stroke:var(--grid)}
:root .figure-host circle{stroke:var(--card)}
:root .graph-stage text{fill:var(--ink)}
:root button:focus-visible{outline:2px solid var(--acc);outline-offset:2px}
:root .annotation-section-tabs{display:flex;gap:18px;margin:0 0 16px;padding:0 4px;border-bottom:1px solid var(--line)}
:root .annotation-section-tabs button{font-size:.8rem!important;font-weight:560;padding:7px 2px 9px!important;border:0;border-bottom:2px solid transparent;border-radius:0;background:transparent!important;color:var(--mut)}
:root .annotation-section-tabs button.active{color:var(--acc)!important;border-bottom-color:var(--acc)}
'''

THEME_SCRIPT = r'''
(function(){
  var root=document.documentElement;
  if(typeof matchMedia!=='function')return;
  var media=matchMedia('(prefers-color-scheme: dark)');
  function update(){
    if(root.hasAttribute('data-theme-manual'))return;
    root.dataset.theme=media.matches?'dark':'light';
    document.dispatchEvent(new CustomEvent('metahopper:theme'));
    if(typeof window.__draw==='function')window.__draw();
  }
  update();media.addEventListener('change',update);
})();
'''
