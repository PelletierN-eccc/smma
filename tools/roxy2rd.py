"""Minimal roxygen-to-Rd converter used to build man/ without R.
(Run devtools::document() instead when R and roxygen2 are available.)"""
import re, glob, os, sys

def md(s):
    s = s.replace('%', '\\%')
    s = re.sub(r'\[([A-Za-z_.][\w.]*)\(\)\]', r'\\code{\\link{\1}()}', s)
    s = re.sub(r'`([^`]+)`', lambda m: '\\code{' + m.group(1).replace('{', '\\{').replace('}', '\\}') + '}', s)
    return s

blocks = []
for f in sorted(glob.glob('R/*.R')):
    lines = open(f).read().split('\n')
    i = 0
    while i < len(lines):
        if lines[i].startswith("#'"):
            j = i
            while j < len(lines) and lines[j].startswith("#'"): j += 1
            doc = [l[2:].lstrip(' ') if l.startswith("#' ") else l[2:] for l in lines[i:j]]
            k = j
            while k < len(lines) and not lines[k].strip(): k += 1
            code = lines[k] if k < len(lines) else ''
            sig = []
            if '<- function' in code:
                m = k
                buf = ''
                while m < len(lines):
                    buf += lines[m].strip() + ' '
                    if buf.count('(') and buf.count('(') == buf.count(')') and '{' in buf: break
                    m += 1
                name = code.split('<-')[0].strip()
                args = buf[buf.index('function') + 8:]
                depth = 0; out = ''
                for ch in args:
                    out += ch
                    if ch == '(': depth += 1
                    if ch == ')':
                        depth -= 1
                        if depth == 0: break
                sig = (name, out.strip())
            elif code.strip().startswith('"_PACKAGE"'):
                sig = ('smma-package', None)
            blocks.append((doc, sig))
            i = j
        else:
            i += 1

parsed = {}
for doc, sig in blocks:
    if not sig: continue
    name, args = sig
    tags = {'param': [], 'return': '', 'examples': '', 'inherit': [], 'export': False, 'keywords': None}
    text, cur = [], None
    for l in doc:
        m = re.match(r'@(\w+)\s*(.*)', l)
        if m:
            t, rest = m.groups()
            cur = t
            if t == 'param':
                pn, pd = (rest.split(' ', 1) + [''])[:2]
                tags['param'].append([pn, pd])
            elif t == 'return': tags['return'] = rest
            elif t == 'examples': tags['examples'] = ''
            elif t == 'inheritParams': tags['inherit'].append(rest.strip())
            elif t == 'export': tags['export'] = True
            elif t == 'keywords': tags['keywords'] = rest.strip()
        else:
            if cur is None: text.append(l)
            elif cur == 'param': tags['param'][-1][1] += '\n' + l
            elif cur == 'return': tags['return'] += '\n' + l
            elif cur == 'examples': tags['examples'] += l + '\n'
    paras = '\n'.join(text).strip().split('\n\n')
    title = paras[0].replace('\n', ' ').strip()
    desc = paras[1:] if len(paras) > 1 else [title]
    parsed[name] = dict(title=title, desc=desc, args=args, **tags)

for name, p in parsed.items():
    params = {}
    for pn, pd in p['param']:
        for n in pn.split(','): params[n.strip()] = (pn, pd)
    for src in p['inherit']:
        for pn, pd in parsed.get(src, {}).get('param', []):
            for n in pn.split(','):
                n = n.strip()
                if p['args'] and re.search(r'[(,]\s*' + re.escape(n) + r'\s*[=,)]', p['args']) and n not in params:
                    params[n] = (n, pd)
    if name == 'smma-package':
        rd = ['\\docType{package}', '\\name{smma-package}', '\\alias{smma}', '\\alias{smma-package}',
              '\\title{' + md(p['title']) + '}', '\\description{', md('\n\n'.join(p['desc'])), '}',
              '\\keyword{internal}']
        open('man/smma-package.Rd', 'w').write('\n'.join(rd) + '\n'); continue
    if not p['export']: continue
    rd = ['% Generated from roxygen comments by tools/roxy2rd.py', '\\name{%s}' % name, '\\alias{%s}' % name,
          '\\title{%s}' % md(p['title'])]
    usage = name + p['args'] if not name.startswith('print.') else '\\method{print}{%s}%s' % (name[6:], p['args'])
    usage = re.sub(r'\s+', ' ', usage)
    rd += ['\\usage{', usage.replace('%', '\\%'), '}']
    if params:
        rd.append('\\arguments{')
        seen = set()
        argnames = re.findall(r'[(,]\s*([A-Za-z_.][\w.]*)\s*(?==|,|\))', p['args'] or '')
        for a in argnames:
            if a in params and params[a][0] not in seen:
                pn, pd = params[a]; seen.add(pn)
                rd.append('\\item{%s}{%s}' % (pn, md(pd.strip())))
            elif a == '...' or a not in params:
                pass
        if '...' in (p['args'] or '') and '...' in params:
            rd.append('\\item{...}{%s}' % md(params['...'][1].strip()))
        rd.append('}')
    rd += ['\\description{', md('\n\n'.join(p['desc'])), '}']
    if p['return']: rd += ['\\value{', md(p['return'].strip()), '}']
    if p['examples'].strip(): rd += ['\\examples{', p['examples'].rstrip().replace('%', '\\%'), '}']
    open('man/%s.Rd' % name.replace('.', '_') if name.startswith('print.') else 'man/%s.Rd' % name, 'w').write('\n'.join(rd) + '\n')
print(len(os.listdir('man')), 'Rd files')
