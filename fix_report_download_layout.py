#!/usr/bin/env python3
"""Patch only annotation download styling in an existing MetaHopper HTML report.

Usage: python fix_report_download_layout.py -i SAMPLE/metahopper_report.html
Edits in place with a one-time .before-download-layout.html backup. No dependencies.
"""
from __future__ import annotations

import argparse
from html.parser import HTMLParser
import os
from pathlib import Path
import re
import stat
import tempfile


STYLE_ID = 'metahopper-download-layout'
DOWNLOAD_LAYOUT_CSS = '''
/* Keep annotation download controls compact in current and older reports. */
#pane-annotation .scroll{max-width:100%;overflow-x:auto}
#pane-annotation table.sortable{width:max-content!important;min-width:100%;table-layout:auto!important}
#pane-annotation table.sortable th,#pane-annotation table.sortable td{padding:6px 10px;line-height:1.35;vertical-align:middle;white-space:nowrap}
#pane-annotation td.annotation-path{width:auto!important;min-width:370px!important;max-width:none!important;white-space:nowrap!important;word-break:normal!important;overflow-wrap:normal!important}
#pane-annotation td.annotation-path .asset-links{display:inline-flex!important;flex-wrap:nowrap!important;gap:6px;align-items:center;min-width:max-content;white-space:nowrap!important}
#pane-annotation td.annotation-path a{display:inline-block;flex:none;max-width:none!important;white-space:nowrap!important;word-break:normal!important;overflow-wrap:normal!important;border:1px solid var(--line,#dae3e7);border-radius:5px;padding:3px 7px;font-size:.79rem;line-height:1.25;color:var(--acc,#136b78);text-decoration:none;background:var(--card,#fff)}
#pane-annotation td.annotation-path a:hover{background:var(--accs,#e2f0f2)}
#pane-annotation td.annotation-path a:focus-visible{outline:2px solid var(--acc,#136b78);outline-offset:2px}
#pane-annotation table.sortable td.annotation-reason,#pane-annotation table.sortable td:nth-last-child(2){white-space:normal;min-width:90px;max-width:260px}
@media(min-width:900px){
 #pane-annotation table.sortable td.annotation-path{position:sticky;right:0;z-index:1;background:var(--card,#fff);box-shadow:-4px 0 7px rgba(22,34,43,.06)}
 #pane-annotation table.sortable th:last-child{position:sticky;right:0;z-index:3;background:var(--accs,#e2f0f2)}
}
'''


class ReportCheck(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.annotation_pane = False
        self.download_cells = 0

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if attrs.get('id') == 'pane-annotation':
            self.annotation_pane = True
        if tag == 'td' and 'annotation-path' in (attrs.get('class') or '').split():
            self.download_cells += 1


def patch_html(text):
    check = ReportCheck()
    check.feed(text)
    check.close()
    if not check.annotation_pane or not check.download_cells:
        raise ValueError('No MetaHopper annotation download table found; no changes made')
    block = '<style id="' + STYLE_ID + '">\n' + DOWNLOAD_LAYOUT_CSS + '</style>'
    existing = re.compile(r'''<style\b(?=[^>]*\bid\s*=\s*(["'])metahopper-download-layout\1)[^>]*>.*?</style\s*>''', re.I | re.S)
    if existing.search(text):
        return existing.sub(lambda _: block, text)
    head = re.search(r'</head\s*>', text, re.I)
    if not head:
        raise ValueError('Report has no closing head element; no changes made')
    return text[:head.start()] + block + '\n' + text[head.start():]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input', type=Path, required=True, help='Existing MetaHopper HTML report')
    parser.add_argument('-o', '--output', type=Path, help='Write a separate patched HTML instead of editing in place')
    args = parser.parse_args()
    temporary = None
    try:
        source = args.input.expanduser().resolve(strict=True)
        original = source.read_bytes()
        patched = patch_html(original.decode('utf-8')).encode('utf-8')
        target = args.output.expanduser().resolve() if args.output else source
        if target == source and patched == original:
            print(f'Layout already fixed: {source}')
            return 0
        if target == source:
            backup = source.with_name(source.stem + '.before-download-layout' + source.suffix)
            try:
                with backup.open('xb') as handle:
                    handle.write(original)
            except FileExistsError:
                pass
            else:
                print(f'Original report backed up: {backup}')
        with tempfile.NamedTemporaryFile(dir=target.parent, prefix='.mh-layout-', delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(patched)
        temporary.chmod(stat.S_IMODE(source.stat().st_mode))
        os.replace(temporary, target)
        temporary = None
        print(f'Fixed download layout: {target}')
        return 0
    except (OSError, ValueError, UnicodeError) as exc:
        parser.exit(1, f'ERROR: {exc}\n')
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


if __name__ == '__main__':
    raise SystemExit(main())
