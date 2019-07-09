#!/usr/bin/env python3

# CharCut: lightweight character-based MT output highlighting and scoring.
# Copyright (C) 2017,2019 Lardilleux
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
CharCut: lightweight character-based MT output highlighting and scoring.

This module can be used as a standalone program, but it also exposes some functionalities
so that they be reused in other projects.
"""

import argparse
import difflib
import gzip
import math
import os.path
import re

from collections import defaultdict, Counter
from itertools import chain
from operator import itemgetter
from xml.sax.saxutils import escape


def make_base_parser():
    """Initiates a CL parser with base arguments used by standalone program and other calling modules"""
    parser = argparse.ArgumentParser(description='''Character-based difference
        highlighting and scoring, based on loose differences. By default, just print the document-level
        score on stdout (0~1, lower is better).''')
    parser.add_argument('-m', '--match-size', type=int, default=3,
                        help='min match size in characters (default: %(default)s)')
    parser.add_argument('-n', '--alt-norm', action='store_true',
                        help='''alternative normalization scheme: use only the candidate's length
                             for normalization (default: %(default)s)''')
    return parser


def add_parser_output_options(parser):
    """Populate a CL base parser with output options"""
    parser.add_argument('-o', '--html-output-file',
                        help='generate a html file with per-segment scores and highlighting')
    parser.add_argument('-p', '--plain-output-file',
                        help='generate a plain text file with per-segment scores only')


def parse_args():
    """Parse and return command line options."""
    parser = make_base_parser()
    add_parser_output_options(parser)
    parser.add_argument('-s', '--src_file', help='source file, only used for display')
    parser.add_argument('file_pair', nargs='+',
                        help='list of comma-separated file pairs to compare, e.g. '
                             'mt1.txt,ref.txt mt2.txt,ref.txt mt1.txt,mt2.txt')
    return parser.parse_args()


def read_gz8(filename):
    """Read a utf8, possibly gzipped, file into memory, as a list of lines."""
    opener = gzip.open if filename.endswith('.gz') else open
    with opener(filename, 'rb') as f:
        return [line.decode('u8') for line in f]


def load_input_files(args):
    """
    Load input files specified in the CL arguments into memory.

    Returns a list of 4-tuples: (segment_id, origin, src_segment,
                                 [(candidate_segment, reference_segment), ...])
    "origin" is always None (present for compatibility with other modules handling sdlxliff files).
    "src_segment" is None if the source file was not passed on the CL.
    There is one (candidate_segment, reference_segment) for each positional argument on the command line.
    """
    cand_ref_file_pairs = [pair.split(',') for pair in args.file_pair]
    cand_ref_zips = [list(map(read_gz8, file_pair)) for file_pair in cand_ref_file_pairs]
    # src file is optional
    src_segs = read_gz8(args.src_file) if args.src_file else [None] * len(cand_ref_zips[0][0])
    for cand_segs, ref_segs in cand_ref_zips:
        assert len(src_segs) == len(cand_segs) == len(ref_segs)
    # Transpose lists
    cand_ref_segs = list(zip(*[list(zip(*cand_ref)) for cand_ref in cand_ref_zips]))
    return [(i, None, src.strip() if src else src, [(cand.strip(), ref.strip()) for cand, ref in cand_refs])
            for i, (src, cand_refs) in enumerate(zip(src_segs, cand_ref_segs), 1)]


def iter_common_substrings(seq1, seq2, start_pos1, start_pos2, min_match_size, add_fix):
    """
    Iterates over common substrings between two sequences, looking at specific start positions.

    start_pos1 (resp. start_pos2) is a list of indices in seq1 (resp. seq2)
    where to look for identical suffixes. This is typically range(len(seq1)) (resp. seq2).

    min_match_size specifies the minimal length of substrings output.

    add_fix is a boolean which indicates whether we should systematically output
    the longest prefix and the longest suffix, independently of min_match_size.

    Returns an iterator over triples: (common substring, [start idx in seq1], [start idx in seq2])
    """
    # We test for equality between elements at positions specified by start_pos1 and start_pos2.
    # For all equal pairs, we perform the test again at matching positions + 1.
    # This is basically a recursive function, but given how inefficient they are in Python,
    # we use a while loop with a stack of "todo" arguments instead.
    # This is simpler than suffix arrays given our constraints (and fast enough).
    n1 = len(seq1)
    n2 = len(seq2)
    # Parameters to the "recursive function". 3rd one is the offset.
    todo = [(start_pos1, start_pos2, 0)]
    while todo:
        pos1, pos2, offset = todo.pop()
        # Associate to each token the list of positions it appears at
        tokens1 = defaultdict(list)
        tokens2 = defaultdict(list)
        for i in pos1:
            if i + offset < n1:
                tokens1[seq1[i+offset]].append(i)
        for i in pos2:
            if i + offset < n2:
                tokens2[seq2[i+offset]].append(i)
        # Take intersection of the two token sets
        for token, ok_pos1 in tokens1.items():
            ok_pos2 = tokens2.get(token)
            if ok_pos2:
                first_pos = ok_pos1[0]
                substr = ''.join(seq1[first_pos:first_pos+offset+1])
                if len(substr) >= min_match_size:
                    yield substr, ok_pos1, ok_pos2
                elif add_fix and 0 in ok_pos1 and 0 in ok_pos2:  # common prefix
                    yield substr, [0], [0]
                elif add_fix and n1-1-offset in ok_pos1 and n2-1-offset in ok_pos2:  # common suffix
                    yield substr, [n1-1-offset], [n2-1-offset]
                todo.append((ok_pos1, ok_pos2, offset+1))


WORD_RE = re.compile('(\W)', re.UNICODE)

def word_split(seq):
    """
    Prepares a sequence of characters for the search of inter-words common substrings.

    The search will be operated directly on the tokens returned (word-based comparison).
    1 non-word character = 1 token.

    Returns an iterator over tuples: (start position, token)

    >>> list(word_split('ab, cd'))
    [(0, 'ab'), (2, ','), (3, ' '), (4, 'cd')]
    """
    pos = 0
    for elt in WORD_RE.split(seq):
        if elt:
            yield pos, elt
            pos += len(elt)


def word_based_matches(seq1, seq2, min_match_size):
    """Iterator over all word-based common substrings between seq1 and seq2."""
    starts1, words1 = list(zip(*word_split(seq1))) if seq1 else ([], [])
    starts2, words2 = list(zip(*word_split(seq2))) if seq2 else ([], [])
    it = iter_common_substrings(words1, words2, list(range(len(words1))),
                                list(range(len(words2))), min_match_size, True)
    for substr, pos1, pos2 in it:
        # Replace positions in words with positions in characters
        yield substr, [starts1[i] for i in pos1], [starts2[i] for i in pos2]


def start_pos(words):
    """Iterator over start positions of a list of words (cumulative lengths)."""
    pos = 0
    for elt in words:
        yield pos
        pos += len(elt)


CHAR_RE = re.compile('(\w+)', re.UNICODE)

def char_split(seq, sep_sign):
    """
    Prepares a sequence of characters for the search of intra-words common substrings.

    Runs of non-word characters are duplicated: they are used once for the preceding word,
    and once for the following word.
    For instance, given the string "ab, cd", we will look for common substrings in
    "ab, " and in ", cd".

    A unique, dummy element is inserted between them in order to prevent the subsequent search
    for common substrings from spanning multiple words. For this purpose, "sep_sign" should be 1
    for the first sequence and -1 for the second one.

    Returns an iterator over triples: (original position, character, is-start-position)
    "is-start-position" is False for trailing non-word characters.

    >>> list(char_split('ab, cd', 1))
    [(0, 'a', True), (1, 'b', True), (2, ',', False), (3, ' ', False),
     (None, 2, False),
     (2, ',', True), (3, ' ', True), (4, 'c', True), (5, 'd', True)]
    """
    split = CHAR_RE.split(seq)
    # Fix in case seq contains only non-word characters
    tokens = ['', split[0], ''] if len(split) == 1 else split
    # "tokens" alternate actual words and runs of non-word characters
    starts = list(start_pos(tokens))
    for i in range(0, len(tokens)-2, 2):
        # insert unique separator to prevent common substrings to span multiple words
        if i:
            yield None, i * sep_sign, False
        for j in range(i, i+3):
            is_start_pos = j != i+2
            for k, char in enumerate(tokens[j], starts[j]):
                yield k, char, is_start_pos


def char_based_matches(seq1, seq2, min_match_size):
    """Iterator over all intra-word character-based common substrings between seq1 and seq2."""
    starts1, chars1, is_start1 = list(zip(*char_split(seq1, 1))) if seq1 else ([], [], [])
    starts2, chars2, is_start2 = list(zip(*char_split(seq2, -1))) if seq2 else ([], [], [])
    start_pos1 = [i for i, is_start in enumerate(is_start1) if is_start]
    start_pos2 = [i for i, is_start in enumerate(is_start2) if is_start]
    ics = iter_common_substrings(chars1, chars2, start_pos1, start_pos2, min_match_size, False)
    for substr, pos1, pos2 in ics:
        # Replace positions with those from the original sequences
        yield substr, [starts1[i] for i in pos1], [starts2[i] for i in pos2]


def order_key(match):
    """Sort key for common substrings: longest first, plus a few heuristic comparisons."""
    substr, pos1, pos2 = match
    return -len(substr), len(pos1) == len(pos2), len(pos1) + len(pos2), pos1


def clean_match_list(match_list, mask1, mask2):
    """
    Filter list of common substrings: remove those for which at least one character
    has already been covered (specified by the two masks).
    """
    for substr, pos1, pos2 in match_list:
        k = len(substr)
        clean_pos1 = [i for i in pos1 if all(mask1[i:i+k])]
        if clean_pos1:
            clean_pos2 = [i for i in pos2 if all(mask2[i:i+k])]
            if clean_pos2:
                yield substr, clean_pos1, clean_pos2


def residual_diff(mask):
    """
    Factor successive 0's from a mask.

    Returns list of pairs: (start position, length)
    """
    buf = []
    for i, elt in enumerate(mask):
        if elt:
            buf.append(i)
        elif buf:
            yield buf[0], len(buf)
            buf = []
    if buf:
        yield buf[0], len(buf)


def greedy_matching(seq1, seq2, min_match_size):
    """
    Greedy search for common substrings between seq1 and seq2.

    Residual substrings (smaller than min_match_size) are also output as deletions (from seq1)
    or insertions (into seq2).

    Returns an iterator over triples: (position in seq1, position in seq2, substring)
    The position in seq1 is -1 for insertions, and the position in seq2 is -1 for deletions.
    """
    assert min_match_size > 0
    retained_matches = []
    # Indicate for each character if it is already covered by a match
    mask1 = [1] * len(seq1)
    mask2 = [1] * len(seq2)

    # List *all* common substrings and sort them (mainly) by length.
    # This is fine since we do (should) not deal with huge strings.
    match_it = chain(word_based_matches(seq1, seq2, min_match_size),
                     char_based_matches(seq1, seq2, min_match_size))
    dedup = {match[0]: match for match in match_it}
    match_list = sorted(dedup.values(), key=order_key)

    # Consume all common substrings, longest first
    while match_list:
        substr, pos1, pos2 = match_list[0]
        i, j = pos1[0], pos2[0]
        retained_matches.append((i, j, substr))
        size = len(substr)
        # Update masks with newly retained characters
        mask1[i:i+size] = [0] * size
        mask2[j:j+size] = [0] * size
        # Eliminate common substrings for which at least one char is already covered
        match_list = list(clean_match_list(match_list, mask1, mask2))

    # Output matches
    for match in retained_matches:
        yield match
    # Output deletions
    for pos, size in residual_diff(mask1):
        yield pos, -1, seq1[pos:pos + size]
    # Output insertions
    for pos, size in residual_diff(mask2):
        yield -1, pos, seq2[pos:pos + size]


def find_regular_matches(ops):
    """
    Find the set of regular (non-shift) matches from the list of operations.

    "ops" is the list of triples as returned by greedy_matching().
    """
    matches1 = sorted(m for m in ops if m[0] != -1 and m[1] != -1)
    matches2 = sorted(matches1, key=lambda match: match[1])
    # Search for the longest common subsequence in characters
    # Expand "string" matches into "character" matches
    char_matches1 = [(m, i) for m in matches1 for i in range(len(m[2]))]
    char_matches2 = [(m, i) for m in matches2 for i in range(len(m[2]))]
    sm = difflib.SequenceMatcher(None, char_matches1, char_matches2, autojunk=False)
    return {m for a, _, size in sm.get_matching_blocks()
            for m, _ in char_matches1[a:a + size]}


def eval_shift_distance(shift, reg_matches):
    """
    Compute the distance in characters a match has been shifted over.

    "reg_matches" is the set of regular matches as returned by find_regular_matches().

    The distance is defined as the number of characters between the shifted match
    and the closest regular match.
    """
    mid_matches = sorted(m for m in reg_matches
                         if (m[0] < shift[0] and m[1] > shift[1])
                         or (m[0] > shift[0] and m[1] < shift[1]))
    return (-(shift[0] - mid_matches[0][0])
            if mid_matches[0][0] < shift[0]
            else (mid_matches[-1][0] + len(mid_matches[-1][2])
                  - (shift[0] + len(shift[2]))))


def add_shift_distance(ops, reg_matches):
    """
    Decorate the list of operations with the shift distance.

    The distance is 0 for everything but shifts.

    Returns an iterator over 4-tuples:
    (pos in seq1, pos in seq2, substring, integer distance)
    """
    # Experimental: turn shifts back into insertions/deletions
    # if the shift distance is "too large".
    for op in ops:
        alo, blo, slice = op
        if alo == -1 or blo == -1 or op in reg_matches:
            yield op + (0,)
        else:  # shift
            dist = eval_shift_distance(op, reg_matches)
            # Heuristic: the shorter a string,
            # the shorter the distance it is allowed to travel
            if math.exp(len(slice)) >= abs(dist):
                yield op + (dist,)
            else:  # replace shift with deletion + insertion
                yield -1, blo, slice, 0
                yield alo, -1, slice, 0


def _merge_adjacent_diffs_aux(diffs):
    prev_start = 0
    prev_substr = ''
    for start, substr in diffs:
        if start == prev_start + len(prev_substr):
            prev_substr += substr
        else:
            if prev_substr:
                yield prev_start, prev_substr
            prev_start = start
            prev_substr = substr
    if prev_substr:
        yield prev_start, prev_substr


def merge_adjacent_diffs(ops):
    """Final cleaning: merge adjacent deletions or insertions into a single operation."""
    matches = [op for op in ops if op[0] != -1 and op[1] != -1]
    deletions = sorted((alo, substr) for alo, blo, substr, _ in ops if blo == -1)
    insertions = sorted((blo, substr) for alo, blo, substr, _ in ops if alo == -1)
    for op in matches:
        yield op
    for alo, substr in _merge_adjacent_diffs_aux(deletions):
        yield alo, -1, substr, 0
    for blo, substr in _merge_adjacent_diffs_aux(insertions):
        yield -1, blo, substr, 0


def add_css_classes(ops):
    """
    Decorate the list of operations with CSS classes for display.

    Each operation is assigned 2 classes:
    * {ins,del,shift,match} for the display style
    * {diff,shift,match}X serve as ids for mouse-overs (substrings that match
    in the two segments compared have the same id)

    Returns an iterator over 6-tuples:
    (pos in seq1, pos in seq2, substring, distance, css class, css id)
    """
    # Substrings are identified based on their start index in the first sequence
    match_alo = 0
    for op in ops:
        alo, blo, _, dist = op
        if alo == -1:
            yield op + ('ins', 'diff{}'.format(match_alo))
        elif blo == -1:
            yield op + ('del', 'diff{}'.format(match_alo))
        elif dist:
            yield op + ('shift', 'shift{}'.format(alo))
        else:
            yield op + ('match', 'match{}'.format(alo))
            match_alo = alo


def compare_segments(cand, ref, min_match_size):
    """
    Main segment comparison function.

    cand and ref are the original unicode strings.

    Return a pair of operation list (same 6-tuples as returned by add_css_classes())
    """
    base_ops = list(greedy_matching(cand, ref, min_match_size))
    reg_matches = find_regular_matches(base_ops)
    clean_ops = list(merge_adjacent_diffs(list(add_shift_distance(base_ops, reg_matches))))
    cand_ops = sorted(op for op in clean_ops if op[0] != -1)
    ref_ops = sorted((op for op in clean_ops if op[1] != -1), key=itemgetter(1))
    styled_cand = list(add_css_classes(cand_ops))
    styled_ref = list(add_css_classes(ref_ops))
    return styled_cand, styled_ref


def _get_cost(styled_ops, css_clazz):
    return sum(len(slice) for _, _, slice, _, clazz, _ in styled_ops
               if clazz == css_clazz)


def score_pair(cand, ref, styled_cand, styled_ref, alt_norm):
    """Score a single candidate/reference pair."""
    ins_cost = _get_cost(styled_cand, 'del')
    del_cost = _get_cost(styled_ref, 'ins')
    # shifts are identical in cand and ref
    shift_cost = _get_cost(styled_cand, 'shift')
    cost = ins_cost + del_cost + shift_cost
    div = 2 * len(cand) if alt_norm else len(cand) + len(ref)
    # Prevent scores > 100%
    bounded_cost = min(cost, div)
    return bounded_cost, div


def score_all(aligned_segs, styled_ops, alt_norm):
    """Score segment pairs based on their differences."""
    for ((_, _, _, cand_refs), styled_cand_refs) in zip(aligned_segs, styled_ops):
        yield [score_pair(cand, ref, styled_cand, styled_ref, alt_norm)
               for (cand, ref), (styled_cand, styled_ref) in zip(cand_refs, styled_cand_refs)]


def ops2html(styled_ops, seg_id):
    for op in styled_ops:
        _, _, slice, dist, css, css_id = op
        substr_id = 'seg{}_{}'.format(seg_id, css_id)
        dist_str = '({:+d})'.format(dist) if dist else ''
        slice_len = len(slice)
        esc_slice = escape(slice)
        yield '<span title="{css}{dist_str}: {slice_len}" class="{css} {substr_id}" ' \
              'onmouseenter="enter(\'{substr_id}\')" onmouseleave="leave(\'{substr_id}\')">' \
              '{esc_slice}</span>'.format(**locals())


def cand2html(ops, seg_id):
    for i, (styled_cand, _) in enumerate(ops):
        cand_id = '{}.{}'.format(seg_id, i)
        yield '''        <td class="compcol trg">{}</td>'''.format(''.join(ops2html(styled_cand, cand_id)))


def ref2html(ops, seg_id):
    for i, (_, styled_ref) in enumerate(ops):
        cand_id = '{}.{}'.format(seg_id, i)
        yield '''        <td class="compcol trg">{}</td>'''.format(''.join(ops2html(styled_ref, cand_id)))


def scores2html(score_pairs):
    for cost, div in score_pairs:
        score = (1. * cost / div) if div else 0
        yield '''        <td class="score"><span class="detail">{:.0f}/{:.0f}= </span>{:.0%}</td>'''.format(cost, div, score)


def basenames2html(idx, basename_pairs, bullets):
    # idx=0 for cand; idx=1 for ref
    cell = '        <td class="compcol seghead">{}{}</td>'
    return '\n'.join(cell.format(bullets.get(pair[idx], ''), pair[idx])
                     for pair in basename_pairs)


def segs2html(segs, ops, score_pairs, src_basename, basename_pairs, bullets):
    """Apply highlighting on a single segment pair."""
    seg_id, origin, src, cand_refs = segs
    origin_str = '<p class="info">({})</p>'.format(origin) if origin else ''
    src_str = '''      <tr><td colspan="{}" class="seghead">{}</td></tr>
      <tr><td colspan="{}" class="srcrow src">{}</td></tr>'''.format(
      len(ops), src_basename, len(ops), escape(src)) if src else '<!-- no source -->'
    cand_basename_cells = basenames2html(0, basename_pairs, bullets)
    ref_basename_cells = basenames2html(1, basename_pairs, bullets)
    cand_cells = '\n'.join(cand2html(ops, seg_id))
    ref_cells = '\n'.join(ref2html(ops, seg_id))
    score_cells = '\n'.join(scores2html(score_pairs))
    return '''
<tr>
  <td class="mainrow">{seg_id}{origin_str}</td>
  <td class="mainrow">
    <table>
{src_str}
      <tr>
{cand_basename_cells}
      </tr>
      <tr>
{cand_cells}
      </tr>
      <tr>
{ref_basename_cells}
      </tr>
      <tr>
{ref_cells}
      </tr>
      <tr>
{score_cells}
      </tr>
    </table>
  </td>
</tr>'''.format(**locals())


def make_coloured_bullets(basenames):
    """
    Create a mapping between filenames and corresponding coloured HTML bullet.
    """
    # Predefined list of colours, used in turn.
    # If the filenames are more numerous, the colours are just re-used.
    colours = ['#00d000', '#d0d000', '#d000d0', '#d0d0d0', '#d00000', '#0000d0', '#00d0d0']
    # Colour only filenames that are used at least twice
    to_colour = [name for name, count in Counter(basenames).most_common() if count > 1]
    bullet = '<span style="color:{}">&#9632;</span> '
    return {to_colour: bullet.format(colours[i % len(colours)])
            for i, to_colour in enumerate(to_colour)}


def html_dump(out_file, aligned_segs, styled_ops, seg_scores, doc_costs, doc_divs, args):
    """
    Apply highlighting on all segments and output them in HTML.

    aligned_segs are the input segments as returned by load_input_files().
    styled_ops are the decorated operations as returned by compare_segments().
    seg_scores are the pairs (cost, div) as returned by score_all().
    """
    html_top = '''<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>charcut output</title>
  <style>
    body {font-family: sans-serif; font-size: 11pt;}
    table, td, th {border-spacing: 0;}
    th {padding: 10px;}
    td {padding: 5px;}
    th {border-top: solid black 2px; font-weight: normal;}
    #key {margin-left: auto; margin-right: auto; text-align: left;}
    .tophead {border-bottom: solid black 1px;}
    .src {font-style: oblique;}
    .trg {font-family: Consolas, monospace;}
    .del {font-weight: bold; color: #f00000;}
    .ins {font-weight: bold; color: #0040ff;}
    .shift {font-weight: bold;}
    .match {}
    .mainrow {border-top: solid black 1px; padding: 1em;}
    .srcrow {padding-bottom: 15px;}
    .compcol {border-left: solid lightgray 2px;}
    .seghead {color: gray; padding-bottom: 0;}
    .score {font-family: Consolas, monospace; font-size: large; text-align: center;}
    .detail {font-size: xx-small; color: gray; text-align: right;}
    .info {font-size: xx-small; color: gray;}
  </style>
  <script>
    function enter(cls) {
      var elts = document.getElementsByClassName(cls);
      for (var i=0; i<elts.length; i++)
        elts[i].style.backgroundColor = "yellow";
    }
    function leave(cls) {
      var elts = document.getElementsByClassName(cls);
      for (var i=0; i<elts.length; i++)
        elts[i].style.backgroundColor = "transparent";
    }
  </script>
</head>
<body>
<table>
<tr>
  <th class="tophead">Seg. id</th>
  <th class="tophead">
    <table id="key">
      <tr>
        <td>Segment comparison:</td>
        <td class="trg">
          <span class="del">Deletion</span><br/>
          <span class="ins">Insertion</span><br/>
          <span class="shift">Shift</span>
        </td>
      </tr>
    </table>
  </th>
</tr>'''
    print(html_top, file=out_file)
    src_basename = os.path.basename(args.src_file) if args.src_file else ''
    basename_pairs = [list(map(os.path.basename, pair.split(','))) for pair in args.file_pair]
    bullets = make_coloured_bullets(chain(*basename_pairs))
    prev_id = None
    for segs, ops, score_pairs in zip(aligned_segs, styled_ops, seg_scores):
        if prev_id:
            # There may be mismatches with unsafe input
            try:
                skipped = int(segs[0]) - int(prev_id) - 1
            except ValueError:
                # Some seg ids contain letters, just ignore
                skipped = None
            if skipped:
                print('''
<tr>
  <td class="info" title="Mismatch - {} seg. skipped">[...]</td>
</tr>'''.format(skipped), file=out_file)

        prev_id = segs[0]
        print(segs2html(segs, ops, score_pairs, src_basename, basename_pairs, bullets), file=out_file)

    score_cell = '<td class="score"><span class="detail">{:.0f}/{:.0f}= </span>{:.0%}</td>'
    score_row = ''.join(score_cell.format(doc_cost, doc_div, (1.*doc_cost/doc_div) if doc_div else 0)
                        for doc_cost, doc_div in zip(doc_costs, doc_divs))
    print('''
<tr>
  <th>Total</th>
  <th>
    <table style="width:100%"><tr>{}</tr></table>
  </th>
</tr>
</table>
</html>'''.format(score_row), file=out_file)


def format_score(cost, div):
    score = (1. * cost / div) if div else 0.
    return '{:.4f} ({}/{})'.format(score, cost, div)


def run_on(aligned_segs, args):
    """
    Main function.

    aligned_seg and args are as returned by load_input_files() and parse_args().
    This way this function can be reused by other modules using different arguments
    or input means.

    Returns the document-level score of the first  (0~1).
    """
    
    styled_ops = [[compare_segments(cand, ref, args.match_size) for cand, ref in cand_refs]
                  for seg_id, _, _, cand_refs in aligned_segs]

    seg_scores = list(score_all(aligned_segs, styled_ops, args.alt_norm))
    pair_scores = list(zip(*seg_scores))

    doc_costs = [sum(cost for cost, _ in pairs) for pairs in pair_scores]
    doc_divs = [sum(div for _, div in pairs) for pairs in pair_scores]

    print('\t'.join(format_score(doc_cost, doc_div)
                    for doc_cost, doc_div in zip(doc_costs, doc_divs)))

    if getattr(args, 'plain_output_file', None):
        with open(args.plain_output_file, 'w') as plain_file:
            for pairs in seg_scores:
                print('\t'.join(format_score(*pair) for pair in pairs), file=plain_file)

    if getattr(args, 'html_output_file', None):
        with open(args.html_output_file, 'w') as html_file:
            html_dump(html_file, aligned_segs, styled_ops, seg_scores, doc_costs, doc_divs, args)

    return (1. * doc_costs[0] / doc_divs[0]) if doc_divs[0] else 0., len(aligned_segs)


if __name__ == '__main__':
    args = parse_args()
    aligned_segs = load_input_files(args)
    run_on(aligned_segs, args)
