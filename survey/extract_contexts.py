import xml.etree.ElementTree as ET, re, json, sys, os
NMA_RE = re.compile(r'network meta[\s\-\u2010\u2011\u2012\u2013]?analys|mixed treatment comparison|multiple treatments? meta[\s\-\u2010\u2011\u2012\u2013]?analys|indirect treatment comparison', re.I)
MA_RE = re.compile(r'meta[\s\-\u2010\u2011\u2012\u2013]?analys', re.I)
COMB_RE = re.compile(r'(updated?|cumulative|new|our own|own) meta[\s\-\u2010\u2011\u2012\u2013]?analys|meta[\s\-\u2010\u2011\u2012\u2013]?analys\w* (including|incorporating|that included|which included) (our|the present|this|the current)|pool(ed|ing) (our|these|the present|the current) (results|data|findings)|(combin|add|incorporat|integrat)\w* (our|these|the present|the current|this trial|the trial)\w*[^.]{0,60}(meta[\s\-\u2010\u2011\u2012\u2013]?analys|previous (trials|studies)|prior (trials|studies)|existing (trials|evidence))|prior (distribution|information)[^.]{0,80}meta[\s\-\u2010\u2011\u2012\u2013]?analys|meta[\s\-\u2010\u2011\u2012\u2013]?analys[^.]{0,80}prior (distribution|information)|(previous|earlier|prior) (trials|studies)[^.]{0,60}(pooled|combined|meta[\s\-\u2010\u2011\u2012\u2013]?analy)', re.I)
def txt(e): return re.sub(r'\s+', ' ', ''.join(e.itertext())).strip() if e is not None else ''
def sents(p): return [s for s in re.split(r'(?<=[.!?])\s+(?=[A-Z\[(])', p) if s]
def parse(path):
    r = ET.parse(path).getroot()
    am = r.find('.//article-meta')
    title = txt(am.find('.//article-title'))
    abst = txt(am.find('.//abstract'))
    # references
    refs = {}; order = []
    for ref in r.iter('ref'):
        rid = ref.get('id'); s = txt(ref)
        order.append(rid); refs[rid] = s
    ma_ids = {k for k, v in refs.items() if MA_RE.search(v)}
    nma_ids = {k for k, v in refs.items() if NMA_RE.search(v)}
    pos = {rid: i for i, rid in enumerate(order)}
    out_ctx = []; comb = []; intro_last = ''
    body = r.find('body')
    if body is None: return None
    # walk sections with titles
    def walk(e, path):
        for ch in e:
            if ch.tag == 'sec':
                ti = txt(ch.find('title')); walk(ch, path + [ti])
            elif ch.tag == 'p':
                ps = txt(ch)
                # cited ids in paragraph incl. ranges
                xs = [x for x in ch.iter('xref') if x.get('ref-type') == 'bibr']
                cited = set()
                for x in xs:
                    for rid in (x.get('rid') or '').split():
                        cited.add(rid)
                # expand ranges: consecutive xrefs separated by dash in text
                raw = ET.tostring(ch, encoding='unicode')
                for m in re.finditer(r'rid="([^"]+)"[^>]*>[^<]*</xref>\s*[–—\-]\s*<xref[^>]*rid="([^"]+)"', raw):
                    a, b = m.group(1).split()[0], m.group(2).split()[0]
                    if a in pos and b in pos and pos[a] < pos[b] and pos[b]-pos[a] < 40:
                        for rid in order[pos[a]:pos[b]+1]: cited.add(rid)
                hit = cited & ma_ids
                sec = ' > '.join(t for t in path if t) or '(no section)'
                if hit:
                    # sentences that cite the MA (approximate: sentences containing an xref label of a hit) -> use paragraph sentences with +-1
                    labels = set()
                    for x in xs:
                        if set((x.get('rid') or '').split()) & hit: labels.add(txt(x))
                    ss = sents(ps); keep = set()
                    for i, s in enumerate(ss):
                        if any(l and l in s for l in labels) or MA_RE.search(s):
                            keep.update({i-1, i, i+1})
                    ctx = ' '.join(ss[i] for i in sorted(keep) if 0 <= i < len(ss))
                    out_ctx.append({'section': sec, 'refs': sorted(hit, key=lambda k: pos.get(k, 0)), 'text': ctx[:1500]})
                for s in sents(ps):
                    if COMB_RE.search(s): comb.append({'section': sec, 'text': s[:500]})
        return
    walk(body, [])
    # last paragraph of the introduction (aims)
    for sec in body.findall('sec'):
        ti = txt(sec.find('title')).lower()
        if 'intro' in ti or 'background' in ti:
            ps = sec.findall('.//p')
            if ps: intro_last = txt(ps[-1])[:1200]
            break
    return {'title': title, 'abstract': abst[:1800],
            'ma_refs': {k: refs[k][:300] for k in sorted(ma_ids, key=lambda k: pos.get(k, 0))},
            'nma_ref_ids': sorted(nma_ids, key=lambda k: pos.get(k, 0)),
            'contexts': out_ctx, 'combine_flags': comb[:10], 'intro_last_paragraph': intro_last}
if __name__ == '__main__':
    ids = json.load(open(sys.argv[1])); stratum = sys.argv[2]; outf = sys.argv[3]
    res = []
    for p in ids:
        f = f'xml/{p}.xml'
        if not os.path.exists(f): res.append({'pmcid': p, 'stratum': stratum, 'error': 'missing'}); continue
        try:
            d = parse(f)
        except Exception as e:
            d = None; err = str(e)
        if d is None: res.append({'pmcid': p, 'stratum': stratum, 'error': 'no body'}); continue
        d.update({'pmcid': p, 'stratum': stratum}); res.append(d)
    json.dump(res, open(outf, 'w'), indent=1, ensure_ascii=False)
    import statistics
    ok = [r for r in res if 'error' not in r]
    print(len(res), 'ok', len(ok), 'mean contexts', statistics.mean(len(r['contexts']) for r in ok), 'no nma ref found', sum(1 for r in ok if stratum=='NMA' and not r['nma_ref_ids']), 'chars', sum(len(json.dumps(r)) for r in res))
