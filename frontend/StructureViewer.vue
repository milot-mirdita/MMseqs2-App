<template>
<div class="structure-panel" v-if="alignments.length > 0 && 'tCa' in alignments[0]">
    <StructureViewerTooltip attach=".structure-panel" />
    <div class="structure-wrapper" ref="structurepanel">
        <table v-if="tmAlignResults" class="tmscore-panel" v-bind="tmPanelBindings">
            <tr>
                <td class="left-cell">TM-Score:</td>
                <td class="right-cell">{{ tmAlignResults.tmScore }}</td>
            </tr>
            <tr>
                <td class="left-cell">RMSD:</td>
                <td class="right-cell">{{ tmAlignResults.rmsd  }}</td>
            </tr>
        </table>
        <StructureViewerToolbar
            :isFullscreen="isFullscreen"
            :isSpinning="isSpinning"
            :showQuery="showQuery"
            :showTarget="showTarget"
            :showArrows="showArrows"
            :disableQueryButton="!hasQuery"
            :disableArrowButton="!hasQuery"
            @makeImage="handleMakeImage"
            @makePDB="handleMakePDB"
            @resetView="handleResetView"
            @toggleFullscreen="handleToggleFullscreen"
            @toggleTarget="handleToggleTarget"
            @toggleQuery="handleToggleQuery"
            @toggleArrows="handleToggleArrows"
            @toggleSpin="handleToggleSpin"
        />
        <div class="structure-viewer" ref="viewport"></div>
    </div>
</div>
</template>

<script>
import StructureViewerTooltip from './StructureViewerTooltip.vue';
import StructureViewerToolbar from './StructureViewerToolbar.vue';
import StructureViewerMixin from './StructureViewerMixin.vue';
import { mockPDB, downloadBlob } from './Utilities.js';
import { pulchra } from 'pulchra-wasm';
import { tmalign, parse as parseTMOutput, parseMatrix as parseTMMatrix } from 'tmalign-wasm';
import { MolScriptBuilder as MS } from 'molstar/lib/mol-script/language/builder';
import { Script } from 'molstar/lib/mol-script/script';
import { parsePDB } from 'molstar/lib/mol-io/reader/pdb/parser';
import { parseCifText } from 'molstar/lib/mol-io/reader/cif/text/parser';
import { StructureSelection } from 'molstar/lib/mol-model/structure';
import { Color } from 'molstar/lib/mol-util/color';

const DEFAULT_QUERY_COLOR = 0x1e88e5;
const DEFAULT_QUERY_UNALIGNED_COLOR = 0xa5cff5;
const DEFAULT_TARGET_COLOR = 0xffc107;
const DEFAULT_TARGET_UNALIGNED_COLOR = 0xffe699;
const DEFAULT_HIGHLIGHT_COLOR = 0xff5252;
const DEFAULT_ARROW_COLOR = 0x9e9e9e;
const MAX_ARROW_PAIRS = 10;
const CHAIN_ID_POOL = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';

// Get chain from structure name like Structure_A
const getChainName = (name) => {
    // HACK FIXME fix AF chain names
    if (/_v[0-9]+$/.test(name)) {
        return 'A';
    }
    let pos = name.lastIndexOf('_');
    if (pos != -1) {
        let match = name.substring(pos + 1);
        return match.length >= 1 ? match[0] : 'A';
    }
    // fallback
    return 'A';
}

const extractChainFromName = (name) => {
    if (!name) return null;
    const trimmed = String(name).trim();
    if (/_v[0-9]+$/.test(trimmed)) {
        return null;
    }
    if (trimmed.length === 1) {
        return trimmed;
    }
    const match = trimmed.match(/(?:_|:|\|)([A-Za-z0-9])$/);
    if (match) {
        return match[1];
    }
    return null;
};

const resolveChainFromName = (name, chainSet, chainList) => {
    const chain = extractChainFromName(name);
    if (!chain) {
        return null;
    }
    if (!chainSet || chainSet.size === 0) {
        return chain;
    }
    if (chainSet.has(chain)) {
        return chain;
    }
    const chainIndex = Number.parseInt(chain, 10);
    if (Number.isFinite(chainIndex)) {
        const zeroIndex = chainIndex - 1;
        if (zeroIndex >= 0 && chainList && zeroIndex < chainList.length) {
            return chainList[zeroIndex];
        }
    }
    return null;
};

const assignAlignmentChains = (alignments, nameGetter, chainSet = null, serialMap = null, lengthKey = null) => {
    const map = new Map();
    const used = new Set();
    const chainList = chainSet && chainSet.size > 0 ? Array.from(chainSet) : null;
    alignments.forEach((alignment, index) => {
        const chain = resolveChainFromName(nameGetter(alignment), chainSet, chainList);
        if (chain) {
            map.set(index, chain);
            used.add(chain);
        }
    });
    if (chainSet && chainSet.size > 0 && serialMap && serialMap.size > 0 && lengthKey) {
        const lengthMap = new Map();
        serialMap.forEach((resnos, chain) => {
            lengthMap.set(chain, resnos.length);
        });
        const candidatesByLength = new Map();
        lengthMap.forEach((len, chain) => {
            if (used.has(chain)) {
                return;
            }
            if (!candidatesByLength.has(len)) {
                candidatesByLength.set(len, []);
            }
            candidatesByLength.get(len).push(chain);
        });
        alignments.forEach((alignment, index) => {
            if (map.has(index)) {
                return;
            }
            const len = Number(alignment[lengthKey]);
            if (!Number.isFinite(len)) {
                return;
            }
            const candidates = candidatesByLength.get(len);
            if (candidates && candidates.length === 1) {
                const chain = candidates[0];
                map.set(index, chain);
                used.add(chain);
                candidatesByLength.delete(len);
            }
        });
    }
    const pool = chainList && chainList.length > 0 ? chainList : CHAIN_ID_POOL.split('');
    let poolIndex = 0;
    alignments.forEach((alignment, index) => {
        if (map.has(index)) {
            return;
        }
        while (poolIndex < pool.length && used.has(pool[poolIndex])) {
            poolIndex += 1;
        }
        const chain = pool[poolIndex] || 'A';
        map.set(index, chain);
        used.add(chain);
        poolIndex += 1;
    });
    return map;
};

const toMolstarColor = (value, fallback) => {
    if (typeof value === 'number' && Number.isFinite(value)) {
        return Color(value);
    }
    if (typeof value === 'string') {
        const trimmed = value.trim();
        if (trimmed.length > 0) {
            const hex = trimmed.startsWith('#') ? trimmed.slice(1) : trimmed;
            const parsed = Number.parseInt(hex, 16);
            if (Number.isFinite(parsed)) {
                return Color(parsed);
            }
        }
    }
    return Color(fallback);
};

const isChainToken = (token) => {
    return typeof token === 'string' && /^[A-Za-z]$/.test(token);
};

const parseAtomLine = (line) => {
    if (!line.startsWith('ATOM') && !line.startsWith('HETATM')) return null;
    let chainId = line.length >= 22 ? line[21] : '';
    let resno = Number.parseInt(line.slice(22, 26).trim(), 10);
    let resname = line.slice(17, 20).trim();
    if (!Number.isFinite(resno) || !resname) {
        const parts = line.trim().split(/\s+/);
        if (parts.length >= 6 && isChainToken(parts[4]) && !Number.isNaN(Number.parseInt(parts[5], 10))) {
            chainId = parts[4];
            resno = Number.parseInt(parts[5], 10);
            resname = parts[3] || resname;
        } else if (parts.length >= 5) {
            resno = Number.parseInt(parts[4], 10);
            resname = parts[3] || resname;
        }
    }
    return { chainId, resno, resname };
};

const buildSerialResidueMap = (pdb) => {
    const map = new Map();
    if (!pdb) {
        return map;
    }
    const lastResidue = new Map();
    for (const line of pdb.split('\n')) {
        const parsed = parseAtomLine(line);
        if (!parsed) {
            continue;
        }
        const chainKey = parsed.chainId && parsed.chainId.trim() ? parsed.chainId.trim() : '_';
        const residueKey = `${chainKey}:${parsed.resno}:${parsed.resname}`;
        if (lastResidue.get(chainKey) === residueKey) {
            continue;
        }
        lastResidue.set(chainKey, residueKey);
        if (!map.has(chainKey)) map.set(chainKey, []);
        map.get(chainKey).push(parsed.resno);
    }
    return map;
};

const buildSerialResidueMapFromPdbTokens = (tokens) => {
    const map = new Map();
    if (!tokens || !tokens.data || !tokens.indices) {
        return map;
    }
    const { data, indices } = tokens;
    const count = tokens.count || 0;
    const lastResidue = new Map();
    for (let i = 0; i < count; i += 1) {
        const start = indices[i * 2];
        const end = indices[i * 2 + 1];
        if (start === undefined || end === undefined) {
            continue;
        }
        const line = data.substring(start, end);
        const parsed = parseAtomLine(line);
        if (!parsed) {
            continue;
        }
        const chainKey = parsed.chainId && parsed.chainId.trim() ? parsed.chainId.trim() : '_';
        const residueKey = `${chainKey}:${parsed.resno}:${parsed.resname}`;
        if (lastResidue.get(chainKey) === residueKey) {
            continue;
        }
        lastResidue.set(chainKey, residueKey);
        if (!map.has(chainKey)) {
            map.set(chainKey, []);
        }
        map.get(chainKey).push(parsed.resno);
    }
    return map;
};

const getChainIdsFromPdbTokens = (tokens) => {
    const chains = new Set();
    if (!tokens || !tokens.data || !tokens.indices) {
        return chains;
    }
    const { data, indices } = tokens;
    const count = tokens.count || 0;
    for (let i = 0; i < count; i += 1) {
        const start = indices[i * 2];
        const end = indices[i * 2 + 1];
        if (start === undefined || end === undefined) {
            continue;
        }
        const line = data.substring(start, end);
        const parsed = parseAtomLine(line);
        if (!parsed) {
            continue;
        }
        if (parsed.chainId && parsed.chainId.trim()) {
            chains.add(parsed.chainId.trim());
        }
    }
    return chains;
};

const extractChainsAndSerialMapFromCif = (cifFile) => {
    const chainIds = new Set();
    const serialMap = new Map();
    const block = cifFile?.blocks?.[0];
    const atomSite = block?.categories?.atom_site || block?.categories?.['atom_site'];
    if (!atomSite) {
        return { chainIds, serialMap };
    }
    const groupField = atomSite.getField('group_PDB');
    const authAsymField = atomSite.getField('auth_asym_id');
    const labelAsymField = atomSite.getField('label_asym_id');
    const authSeqField = atomSite.getField('auth_seq_id');
    const labelSeqField = atomSite.getField('label_seq_id');
    const authCompField = atomSite.getField('auth_comp_id');
    const labelCompField = atomSite.getField('label_comp_id');
    const lastResidue = new Map();
    const rowCount = atomSite.rowCount || 0;
    for (let i = 0; i < rowCount; i += 1) {
        const group = (groupField ? groupField.str(i) : '') || 'ATOM';
        if (group !== 'ATOM' && group !== 'HETATM') {
            continue;
        }
        const chainIdRaw = (authAsymField?.str(i) || labelAsymField?.str(i) || '').trim();
        const chainKey = chainIdRaw || '_';
        const resnoStr = (authSeqField?.str(i) || labelSeqField?.str(i) || '').trim();
        const resno = Number.parseInt(resnoStr, 10);
        if (!Number.isFinite(resno)) {
            continue;
        }
        const resname = (authCompField?.str(i) || labelCompField?.str(i) || 'UNK').trim();
        const residueKey = `${chainKey}:${resno}:${resname}`;
        if (lastResidue.get(chainKey) === residueKey) {
            continue;
        }
        lastResidue.set(chainKey, residueKey);
        if (!serialMap.has(chainKey)) {
            serialMap.set(chainKey, []);
        }
        serialMap.get(chainKey).push(resno);
        if (chainIdRaw) {
            chainIds.add(chainIdRaw);
        }
    }
    return { chainIds, serialMap };
};

const readQueryDataWithMolIo = async (queryPdb, format) => {
    if (format === 'pdb') {
        const parsed = await parsePDB(queryPdb, 'query').run();
        if (parsed.isError) {
            return null;
        }
        const tokens = parsed.result?.lines;
        return {
            chainIds: getChainIdsFromPdbTokens(tokens),
            serialMap: buildSerialResidueMapFromPdbTokens(tokens),
        };
    }
    if (format === 'mmcif') {
        const parsed = await parseCifText(queryPdb).run();
        if (parsed.isError) {
            return null;
        }
        return extractChainsAndSerialMapFromCif(parsed.result);
    }
    return null;
};

const mapRangesToAuth = (ranges, serialMap) => {
    if (!serialMap || serialMap.length === 0) {
        return ranges;
    }
    return ranges.map(range => {
        const startIdx = Math.max(1, Math.round(range.start));
        const endIdx = Math.max(1, Math.round(range.end));
        const startPos = Math.min(serialMap.length, startIdx) - 1;
        const endPos = Math.min(serialMap.length, endIdx) - 1;
        const startResno = serialMap[startPos];
        const endResno = serialMap[endPos];
        if (!Number.isFinite(startResno) || !Number.isFinite(endResno)) return null;
        return {
            start: Math.min(startResno, endResno),
            end: Math.max(startResno, endResno),
        };
    }).filter(Boolean);
};

const mapSerialPosition = (serialPos, serialMap) => {
    if (!Number.isFinite(serialPos)) {
        return null;
    }
    if (!serialMap || serialMap.length === 0) {
        return serialPos;
    }
    const index = Math.max(1, Math.round(serialPos));
    const pos = Math.min(serialMap.length, index) - 1;
    const resno = serialMap[pos];
    return Number.isFinite(resno) ? resno : serialPos;
};

const buildChainExpression = (chain, ranges, useChainTest) => {
    const groupBy = MS.struct.atomProperty.macromolecular.residueKey();
    const chainTest = useChainTest
        ? MS.core.rel.eq([
            MS.struct.atomProperty.macromolecular.auth_asym_id(),
            chain,
        ])
        : null;
    if (!ranges || ranges.length === 0) {
        return chainTest
            ? MS.struct.generator.atomGroups({ 'chain-test': chainTest, 'group-by': groupBy })
            : MS.struct.generator.all();
    }
    const rangeExpressions = ranges.map(range => MS.struct.generator.atomGroups({
        ...(chainTest ? { 'chain-test': chainTest } : {}),
        'residue-test': MS.core.rel.inRange([
            MS.struct.atomProperty.macromolecular.auth_seq_id(),
            range.start,
            range.end,
        ]),
        'group-by': groupBy,
    }));
    return rangeExpressions.length === 1
        ? rangeExpressions[0]
        : MS.struct.combinator.merge(rangeExpressions);
};

const getSelectionLoci = (expression, structureRef) => {
    const data = structureRef?.cell?.obj?.data;
    if (!data || !expression) {
        return null;
    }
    const selection = Script.getStructureSelection(expression, data);
    return StructureSelection.toLociWithSourceUnits(selection);
};

const getSelectionLociSafe = (expression, structureRef) => {
    const data = structureRef?.cell?.obj?.data;
    if (!data || !expression) {
        return null;
    }
    const selection = Script.getStructureSelection(expression, data);
    if (StructureSelection.isEmpty(selection)) {
        return null;
    }
    return StructureSelection.toLociWithSourceUnits(selection);
};

const buildResidueAtomExpression = (chain, resno, useChainTest, atomId = 'CA') => {
    const groupBy = MS.struct.atomProperty.macromolecular.residueKey();
    const params = {
        'residue-test': MS.core.rel.eq([
            MS.struct.atomProperty.macromolecular.auth_seq_id(),
            resno,
        ]),
        'group-by': groupBy,
    };
    if (useChainTest) {
        params['chain-test'] = MS.core.rel.eq([
            MS.struct.atomProperty.macromolecular.auth_asym_id(),
            chain,
        ]);
    }
    if (atomId) {
        params['atom-test'] = MS.core.rel.eq([
            MS.struct.atomProperty.macromolecular.label_atom_id(),
            atomId,
        ]);
    }
    return MS.struct.generator.atomGroups(params);
};

const extractAtomLines = (pdb) => {
    if (!pdb) {
        return [];
    }
    return pdb.split('\n').filter(line => line.startsWith('ATOM'));
};

const mergePdbChunks = (chunks) => {
    let serial = 1;
    const lines = [];
    chunks.forEach(chunk => {
        const atomLines = extractAtomLines(chunk);
        atomLines.forEach(line => {
            const padded = line.padEnd(80, ' ');
            const serialStr = serial.toString().padStart(5);
            lines.push(`${padded.slice(0, 6)}${serialStr}${padded.slice(11)}`);
            serial++;
        });
    });
    return lines.length > 0 ? `${lines.join('\n')}\nEND\n` : '';
};

const transformPdb = (pdb, t, u) => {
    if (!pdb) {
        return pdb;
    }
    return pdb.split('\n').map(line => {
        if (!line.startsWith('ATOM') && !line.startsWith('HETATM')) {
            return line;
        }
        const x = Number.parseFloat(line.slice(30, 38));
        const y = Number.parseFloat(line.slice(38, 46));
        const z = Number.parseFloat(line.slice(46, 54));
        if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)) {
            return line;
        }
        const nx = t[0] + u[0][0] * x + u[0][1] * y + u[0][2] * z;
        const ny = t[1] + u[1][0] * x + u[1][1] * y + u[1][2] * z;
        const nz = t[2] + u[2][0] * x + u[2][1] * y + u[2][2] * z;
        const prefix = line.slice(0, 30);
        const suffix = line.slice(54);
        return `${prefix}${nx.toFixed(3).padStart(8)}${ny.toFixed(3).padStart(8)}${nz.toFixed(3).padStart(8)}${suffix}`;
    }).join('\n');
};

const applyChainId = (pdb, chainId) => {
    if (!pdb) {
        return pdb;
    }
    const chain = (chainId && chainId.length > 0) ? chainId[0] : 'A';
    return pdb.split('\n').map(line => {
        if (!line.startsWith('ATOM') && !line.startsWith('HETATM')) {
            return line;
        }
        const padded = line.padEnd(80, ' ');
        return `${padded.slice(0, 21)}${chain}${padded.slice(22)}`;
    }).join('\n');
};

const getChainIdsFromPdb = (pdb) => {
    const chains = new Set();
    if (!pdb) {
        return chains;
    }
    const lines = pdb.split('\n');
    for (const line of lines) {
        const parsed = parseAtomLine(line);
        if (!parsed) continue;
        if (parsed.chainId && parsed.chainId.trim()) {
            chains.add(parsed.chainId.trim());
        }
    }
    return chains;
};

const buildAlignmentRanges = (alignments, kind, chainResolver) => {
    const rangesByChain = new Map();
    alignments.forEach((alignment, index) => {
        const chain = chainResolver
            ? chainResolver(alignment, index)
            : getChainName(kind === 'query' ? alignment.query : alignment.target);
        const start = Number(kind === 'query' ? alignment.qStartPos : alignment.dbStartPos);
        const end = Number(kind === 'query' ? alignment.qEndPos : alignment.dbEndPos);
        if (!Number.isFinite(start) || !Number.isFinite(end)) return;
        if (!rangesByChain.has(chain)) {
            rangesByChain.set(chain, []);
        }
        rangesByChain.get(chain).push({ start, end });
    });
    return rangesByChain;
};

const buildAlignedPairsFromBacktrace = (alignment, qChain, tChain) => {
    const pairs = [];
    const backtrace = alignment.backtrace || '';
    if (!backtrace) {
        return pairs;
    }
    let qPos = Number(alignment.qStartPos);
    let tPos = Number(alignment.dbStartPos);
    if (!Number.isFinite(qPos) || !Number.isFinite(tPos)) {
        return pairs;
    }
    for (let i = 0; i < backtrace.length; i += 1) {
        const code = backtrace[i];
        if (code === 'M') {
            pairs.push({ qChain, tChain, qPos, tPos });
            qPos += 1;
            tPos += 1;
            continue;
        }
        if (code === 'I') {
            qPos += 1;
            continue;
        }
        if (code === 'D') {
            tPos += 1;
        }
    }
    return pairs;
};

const makeSubPdbFromRanges = (pdb, rangesByChain) => {
    if (!pdb) {
        return '';
    }
    const lines = pdb.split('\n');
    const atomLines = lines.filter(line => line.startsWith('ATOM') || line.startsWith('HETATM'));
    const selected = [];
    const allRanges = [];
    rangesByChain.forEach(ranges => {
        ranges.forEach(range => allRanges.push(range));
    });
    const residueCounters = new Map();
    const lastResidue = new Map();
    for (const line of lines) {
        const parsed = parseAtomLine(line);
        if (!parsed) continue;
        const chainKey = parsed.chainId && parsed.chainId.trim() ? parsed.chainId.trim() : '_';
        const residueKey = `${chainKey}:${parsed.resno}:${parsed.resname}`;
        let serialIndex = residueCounters.get(chainKey) || 0;
        if (lastResidue.get(chainKey) !== residueKey) {
            serialIndex += 1;
            residueCounters.set(chainKey, serialIndex);
            lastResidue.set(chainKey, residueKey);
        }
        const ranges = rangesByChain.get(parsed.chainId) || allRanges;
        if (!ranges || ranges.length === 0) continue;
        for (const range of ranges) {
            if (serialIndex >= range.start && serialIndex <= range.end) {
                selected.push(line);
                break;
            }
        }
    }
    if (selected.length === 0) {
        return atomLines.join('\n');
    }
    return selected.join('\n');
};

export default {
    name: "StructureViewer",
    components: {
        StructureViewerTooltip,
        StructureViewerToolbar,
    },
    mixins: [
        StructureViewerMixin,
    ],
    data() {
        return {
            selection: null,
            showArrows: false,
            showQuery: 0,
            showTarget: 0,
            tmAlignResults: null,
            hasQuery: true,
            queryData: null,
            targetData: null,
            lastQueryPdb: '',
            lastTargetPdb: '',
            queryChains: new Set(),
            targetChains: new Set(),
            querySerialMap: new Map(),
            targetSerialMap: new Map(),
            queryAlignmentChains: new Map(),
            targetAlignmentChains: new Map(),
            queryStructureRef: null,
            targetStructureRef: null,
            queryVisualComponentRefs: [],
            targetVisualComponentRefs: [],
            highlightComponentRef: null,
            arrowRefs: [],
            renderToken: 0,
        }
    },
    props: {
        alignments: { type: Array },
        highlights: { type: Array },
        queryFile: { type: String },
        queryAlignedColor: { type: String, default: "#1E88E5" },
        queryUnalignedColor: { type: String, default: "#A5CFF5" },
        targetAlignedColor: { type: String, default: "#FFC107" },
        targetUnalignedColor: { type: String, default: "#FFE699" },
        qRepr: { type: String, default: "cartoon" },
        tRepr: { type: String, default: "cartoon" },
        hits: { type: Object },
        autoViewTime: { type: Number, default: 100 }
    },
    methods: {
        updateAlignmentChains() {
            if (!this.alignments || this.alignments.length === 0) {
                this.queryAlignmentChains = new Map();
                this.targetAlignmentChains = new Map();
                return;
            }
            this.queryAlignmentChains = assignAlignmentChains(this.alignments, (alignment) => alignment.query);
            this.targetAlignmentChains = assignAlignmentChains(this.alignments, (alignment) => alignment.target);
        },
        updateQueryAlignmentChainsFromStructure() {
            if (!this.alignments || this.alignments.length === 0) {
                this.queryAlignmentChains = new Map();
                return;
            }
            if (!this.queryChains || this.queryChains.size === 0) {
                return;
            }
            this.queryAlignmentChains = assignAlignmentChains(
                this.alignments,
                (alignment) => alignment.query,
                this.queryChains,
                this.querySerialMap,
                'qLen',
            );
        },
        updateTargetAlignmentChainsFromStructure() {
            if (!this.alignments || this.alignments.length === 0) {
                this.targetAlignmentChains = new Map();
                return;
            }
            if (!this.targetChains || this.targetChains.size === 0) {
                return;
            }
            this.targetAlignmentChains = assignAlignmentChains(
                this.alignments,
                (alignment) => alignment.target,
                this.targetChains,
                this.targetSerialMap,
                'dbLen',
            );
        },
        getAlignmentChain(kind, alignment, index) {
            const map = kind === 'query' ? this.queryAlignmentChains : this.targetAlignmentChains;
            if (map && map.has(index)) {
                return map.get(index);
            }
            return getChainName(kind === 'query' ? alignment.query : alignment.target);
        },
        handleToggleArrows() {
            if (!this.stage) {
                return;
            }
            this.showArrows = !this.showArrows;
        },
        handleToggleQuery() {
            if (!this.stage) {
                return;
            }
            if (__LOCAL__) {
                this.showQuery = (this.showQuery === 0) ? 1 : 0;
            } else {
                this.showQuery = (this.showQuery === 2) ? 0 : this.showQuery + 1;
            }
        },
        handleResetView() {
            if (!this.stage) {
                return;
            }
            this.focusSelection();
        },
        handleToggleTarget() {
            if (!this.stage) {
                return;
            }
            if (__LOCAL__) {
                this.showTarget = (this.showTarget === 0) ? 1 : 0;
            } else {
                this.showTarget = (this.showTarget === 2) ? 0 : this.showTarget + 1; 
            }
        },
        async handleMakeImage() {
            if (!this.stage) {
                return;
            }
            const wasSpinning = this.isSpinning;
            this.isSpinning = false;
            const title = this.alignments.map(aln => this.hasQuery ? `${aln.query}-${aln.target}` : aln.target).join("_");
            const blob = await this.stage.makeImage();
            if (blob) {
                downloadBlob(blob, `${title}.png`);
            }
            this.isSpinning = wasSpinning;
        },
        handleMakePDB() {
            const qPDB = extractAtomLines(this.lastQueryPdb);
            const tPDB = extractAtomLines(this.lastTargetPdb);
            if (qPDB.length === 0 && tPDB.length === 0) return;
            const title = this.alignments.map(aln => qPDB.length > 0 ? `${aln.query}-${aln.target}` : aln.target);
            let result = null;
            if (qPDB.length > 0 && tPDB.length > 0) {
                result =
`TITLE     ${title.join(" ")}
REMARK     This file was generated by the Foldseek webserver:
REMARK       https://search.foldseek.com
REMARK     Please cite:
REMARK       https://doi.org/10.1101/2022.02.07.479398
MODEL        1
${qPDB.join('\n')}
ENDMDL
MODEL        2
${tPDB.join('\n')}
ENDMDL
END
`
            } else {
                result =
`TITLE     ${title.join(" ")}
REMARK     This file was generated by the Foldseek webserver:
REMARK       https://search.foldseek.com
REMARK     Please cite:
REMARK       https://doi.org/10.1101/2022.02.07.479398
MODEL        1
${tPDB.join('\n')}
ENDMDL
END
`
            }
            downloadBlob(new Blob([result], { type: 'text/plain' }), (title.join("_") + ".pdb"));
        },
        buildSelectionExpression(kind, mode) {
            if (!this.alignments || this.alignments.length === 0) {
                return MS.struct.generator.all();
            }
            if (mode === 2) {
                return MS.struct.generator.all();
            }
            const useRanges = mode === 0;
            const chainSet = kind === 'query' ? this.queryChains : this.targetChains;
            const serialMap = kind === 'query' ? this.querySerialMap : this.targetSerialMap;
            const rangesByChain = new Map();
            this.alignments.forEach((alignment, index) => {
                const chainName = this.getAlignmentChain(kind, alignment, index);
                if (!rangesByChain.has(chainName)) {
                    rangesByChain.set(chainName, []);
                }
                if (useRanges) {
                    const start = Number(kind === 'query' ? alignment.qStartPos : alignment.dbStartPos);
                    const end = Number(kind === 'query' ? alignment.qEndPos : alignment.dbEndPos);
                    if (Number.isFinite(start) && Number.isFinite(end)) {
                        rangesByChain.get(chainName).push({ start, end });
                    }
                }
            });
            const chainExpressions = [];
            rangesByChain.forEach((ranges, chain) => {
                let mappedRanges = ranges;
                if (useRanges) {
                    let chainMap = serialMap.get(chain);
                    if (!chainMap && serialMap.size === 1) {
                        chainMap = Array.from(serialMap.values())[0];
                    }
                    mappedRanges = mapRangesToAuth(ranges, chainMap);
                }
                const useChain = chainSet && chainSet.size > 0 && chainSet.has(chain);
                chainExpressions.push(buildChainExpression(chain, useRanges ? mappedRanges : null, useChain));
            });
            if (chainExpressions.length === 0) {
                return MS.struct.generator.all();
            }
            return chainExpressions.length === 1
                ? chainExpressions[0]
                : MS.struct.combinator.merge(chainExpressions);
        },
        buildAlignedPairs(maxPairs = MAX_ARROW_PAIRS) {
            const pairs = [];
            if (!this.alignments || this.alignments.length === 0) {
                return pairs;
            }
            this.alignments.forEach((alignment, index) => {
                const qChain = this.getAlignmentChain('query', alignment, index);
                const tChain = this.getAlignmentChain('target', alignment, index);
                const qAln = alignment.qAln || '';
                const tAln = alignment.dbAln || alignment.tAln || '';
                if (!qAln || !tAln) {
                    pairs.push(...buildAlignedPairsFromBacktrace(alignment, qChain, tChain));
                    return;
                }
                let qPos = Number(alignment.qStartPos);
                let tPos = Number(alignment.dbStartPos);
                if (!Number.isFinite(qPos) || !Number.isFinite(tPos)) {
                    return;
                }
                const len = Math.min(qAln.length, tAln.length);
                for (let i = 0; i < len; i += 1) {
                    const qChar = qAln[i];
                    const tChar = tAln[i];
                    const qGap = qChar === '-';
                    const tGap = tChar === '-';
                    if (!qGap && !tGap) {
                        pairs.push({ qChain, tChain, qPos, tPos });
                    }
                    if (!qGap) {
                        qPos += 1;
                    }
                    if (!tGap) {
                        tPos += 1;
                    }
                }
            });
            if (!maxPairs || pairs.length <= maxPairs) {
                return pairs;
            }
            const stride = Math.ceil(pairs.length / maxPairs);
            const sampled = [];
            for (let i = 0; i < pairs.length; i += stride) {
                sampled.push(pairs[i]);
            }
            return sampled;
        },
        getMappedResno(kind, chain, serialPos) {
            const serialMap = kind === 'query' ? this.querySerialMap : this.targetSerialMap;
            let chainMap = serialMap.get(chain);
            if (!chainMap && serialMap.size === 1) {
                chainMap = Array.from(serialMap.values())[0];
            }
            return mapSerialPosition(serialPos, chainMap);
        },
        getResidueLoci(structureRef, chain, resno, chainSet) {
            if (!structureRef || !Number.isFinite(resno)) {
                return null;
            }
            const useChain = chainSet && chainSet.size > 0 && chainSet.has(chain);
            const expr = buildResidueAtomExpression(chain, resno, useChain, 'CA');
            return getSelectionLociSafe(expr, structureRef);
        },
        async removeStageRef(ref) {
            if (!this.stage || !ref) return;
            await this.stage.remove(ref);
        },
        async clearStageRefs(refs) {
            if (!Array.isArray(refs) || refs.length === 0) return;
            for (const ref of refs) {
                await this.removeStageRef(ref);
            }
            refs.splice(0, refs.length);
        },
        async clearArrowMeasurements() {
            await this.clearStageRefs(this.arrowRefs);
        },
        async clearQueryVisualComponents() {
            await this.clearStageRefs(this.queryVisualComponentRefs);
        },
        async clearTargetVisualComponents() {
            await this.clearStageRefs(this.targetVisualComponentRefs);
        },
        async clearVisualComponents() {
            await this.clearTargetVisualComponents();
            await this.clearQueryVisualComponents();
            await this.removeStageRef(this.highlightComponentRef);
            this.highlightComponentRef = null;
            await this.clearArrowMeasurements();
        },
        async rerenderArrows() {
            await this.clearArrowMeasurements();
            if (!this.stage || !this.showArrows) {
                return;
            }
            if (!this.queryStructureRef || !this.targetStructureRef) {
                return;
            }
            const pairs = this.buildAlignedPairs(MAX_ARROW_PAIRS);
            if (pairs.length === 0) {
                return;
            }
            for (const pair of pairs) {
                const qResno = this.getMappedResno('query', pair.qChain, pair.qPos);
                const tResno = this.getMappedResno('target', pair.tChain, pair.tPos);
                let qLoci = this.getResidueLoci(this.queryStructureRef, pair.qChain, qResno, this.queryChains);
                let tLoci = this.getResidueLoci(this.targetStructureRef, pair.tChain, tResno, this.targetChains);
                if (!qLoci) {
                    qLoci = getSelectionLociSafe(buildResidueAtomExpression(pair.qChain, qResno, false, 'CA'), this.queryStructureRef)
                        || getSelectionLociSafe(buildResidueAtomExpression(pair.qChain, qResno, false, null), this.queryStructureRef);
                }
                if (!tLoci) {
                    tLoci = getSelectionLociSafe(buildResidueAtomExpression(pair.tChain, tResno, false, 'CA'), this.targetStructureRef)
                        || getSelectionLociSafe(buildResidueAtomExpression(pair.tChain, tResno, false, null), this.targetStructureRef);
                }
                if (!qLoci || !tLoci) {
                    continue;
                }
                const ref = await this.stage.addDistance(qLoci, tLoci);
                if (ref) {
                    this.arrowRefs.push(ref);
                }
            }
        },
        async rerenderHighlight(highlightColor = toMolstarColor(DEFAULT_HIGHLIGHT_COLOR, DEFAULT_HIGHLIGHT_COLOR)) {
            await this.removeStageRef(this.highlightComponentRef);
            this.highlightComponentRef = null;
            if (!this.targetStructureRef) {
                return;
            }
            const highlightExpr = this.buildHighlightExpression();
            if (!highlightExpr) {
                return;
            }
            const highlightComponent = await this.stage.createComponentFromExpression(
                this.targetStructureRef,
                highlightExpr,
                `target-highlight-${this.renderToken}`,
            );
            if (!highlightComponent) {
                return;
            }
            await this.stage.addRepresentation(highlightComponent, {
                type: 'ball-and-stick',
                color: 'uniform',
                colorParams: { value: highlightColor },
            });
            this.highlightComponentRef = highlightComponent;
        },
        isMultimerMode() {
            const first = this.alignments?.[0];
            if (!first || typeof first !== 'object') {
                return false;
            }
            return Object.prototype.hasOwnProperty.call(first, 'complexu')
                && Object.prototype.hasOwnProperty.call(first, 'complext');
        },
        async addSurfaceRepresentation(component, color, alpha = 0.12) {
            await this.stage.addRepresentation(
                component,
                {
                    type: 'molecular-surface',
                    color: 'uniform',
                    colorParams: { value: color },
                },
                { initialState: { alphaFactor: alpha } },
            );
        },
        async renderTargetVisuals(token, targetAlignedColor, targetUnalignedColor, isMultimer) {
            await this.clearTargetVisualComponents();
            if (!this.targetStructureRef) {
                return;
            }
            const targetStructure = this.targetStructureRef;
            const alignedExpr = this.buildSelectionExpression('target', 0);
            const chainExpr = this.buildSelectionExpression('target', 1);
            const fullExpr = MS.struct.generator.all();

            if (this.showTarget === 0) {
                const alignedComponent = await this.stage.createComponentFromExpression(targetStructure, alignedExpr, `target-aligned-${token}`);
                if (alignedComponent) {
                    await this.stage.addRepresentation(alignedComponent, {
                        type: 'cartoon',
                        color: 'uniform',
                        colorParams: { value: targetAlignedColor },
                    });
                    this.targetVisualComponentRefs.push(alignedComponent);
                }
                if (isMultimer) {
                    const surfaceComponent = await this.stage.createComponentFromExpression(targetStructure, alignedExpr, `target-surface-aligned-${token}`);
                    if (surfaceComponent) {
                        await this.addSurfaceRepresentation(surfaceComponent, targetAlignedColor, 0.14);
                        this.targetVisualComponentRefs.push(surfaceComponent);
                    }
                }
                return;
            }

            const visibleExpr = this.showTarget === 2 ? fullExpr : chainExpr;
            const unalignedExpr = MS.struct.modifier.exceptBy({ 0: visibleExpr, by: alignedExpr });
            const unalignedComponent = await this.stage.createComponentFromExpression(targetStructure, unalignedExpr, `target-unaligned-${token}`);
            if (unalignedComponent) {
                await this.stage.addRepresentation(unalignedComponent, {
                    type: 'cartoon',
                    color: 'uniform',
                    colorParams: { value: targetUnalignedColor },
                });
                this.targetVisualComponentRefs.push(unalignedComponent);
            }
            const alignedComponent = await this.stage.createComponentFromExpression(targetStructure, alignedExpr, `target-aligned-${token}`);
            if (alignedComponent) {
                await this.stage.addRepresentation(alignedComponent, {
                    type: 'cartoon',
                    color: 'uniform',
                    colorParams: { value: targetAlignedColor },
                });
                this.targetVisualComponentRefs.push(alignedComponent);
            }
            if (isMultimer) {
                const surfaceSourceExpr = this.showTarget === 2 ? fullExpr : chainExpr;
                const surfaceComponent = await this.stage.createComponentFromExpression(targetStructure, surfaceSourceExpr, `target-surface-${token}`);
                if (surfaceComponent) {
                    await this.addSurfaceRepresentation(surfaceComponent, targetAlignedColor, this.showTarget === 2 ? 0.10 : 0.14);
                    this.targetVisualComponentRefs.push(surfaceComponent);
                }
            }
        },
        async renderQueryVisuals(token, queryAlignedColor, queryUnalignedColor, isMultimer) {
            await this.clearQueryVisualComponents();
            if (!this.queryStructureRef) {
                return;
            }
            const queryStructure = this.queryStructureRef;
            const alignedExpr = this.buildSelectionExpression('query', 0);
            const chainExpr = this.buildSelectionExpression('query', 1);
            const fullExpr = MS.struct.generator.all();

            if (this.showQuery === 0) {
                const alignedComponent = await this.stage.createComponentFromExpression(queryStructure, alignedExpr, `query-aligned-${token}`);
                if (alignedComponent) {
                    await this.stage.addRepresentation(alignedComponent, {
                        type: 'cartoon',
                        color: 'uniform',
                        colorParams: { value: queryAlignedColor },
                    });
                    this.queryVisualComponentRefs.push(alignedComponent);
                }
                if (isMultimer) {
                    const surfaceComponent = await this.stage.createComponentFromExpression(queryStructure, alignedExpr, `query-surface-aligned-${token}`);
                    if (surfaceComponent) {
                        await this.addSurfaceRepresentation(surfaceComponent, queryAlignedColor, 0.14);
                        this.queryVisualComponentRefs.push(surfaceComponent);
                    }
                }
                return;
            }

            const visibleExpr = this.showQuery === 2 ? fullExpr : chainExpr;
            const unalignedExpr = MS.struct.modifier.exceptBy({ 0: visibleExpr, by: alignedExpr });
            const unalignedComponent = await this.stage.createComponentFromExpression(queryStructure, unalignedExpr, `query-unaligned-${token}`);
            if (unalignedComponent) {
                await this.stage.addRepresentation(unalignedComponent, {
                    type: 'cartoon',
                    color: 'uniform',
                    colorParams: { value: queryUnalignedColor },
                });
                this.queryVisualComponentRefs.push(unalignedComponent);
            }
            const alignedComponent = await this.stage.createComponentFromExpression(queryStructure, alignedExpr, `query-aligned-${token}`);
            if (alignedComponent) {
                await this.stage.addRepresentation(alignedComponent, {
                    type: 'cartoon',
                    color: 'uniform',
                    colorParams: { value: queryAlignedColor },
                });
                this.queryVisualComponentRefs.push(alignedComponent);
            }
            if (isMultimer) {
                const surfaceSourceExpr = this.showQuery === 2 ? fullExpr : chainExpr;
                const surfaceComponent = await this.stage.createComponentFromExpression(queryStructure, surfaceSourceExpr, `query-surface-${token}`);
                if (surfaceComponent) {
                    await this.addSurfaceRepresentation(surfaceComponent, queryAlignedColor, this.showQuery === 2 ? 0.10 : 0.14);
                    this.queryVisualComponentRefs.push(surfaceComponent);
                }
            }
        },
        async rerenderTargetVisuals() {
            if (!this.stage || !this.stageReady) {
                return;
            }
            const token = ++this.renderToken;
            await this.stageReady;
            if (token !== this.renderToken) {
                return;
            }
            const targetAlignedColor = toMolstarColor(this.targetAlignedColor, DEFAULT_TARGET_COLOR);
            const targetUnalignedColor = toMolstarColor(this.targetUnalignedColor, DEFAULT_TARGET_UNALIGNED_COLOR);
            await this.renderTargetVisuals(token, targetAlignedColor, targetUnalignedColor, this.isMultimerMode());
        },
        async rerenderQueryVisuals() {
            if (!this.stage || !this.stageReady) {
                return;
            }
            const token = ++this.renderToken;
            await this.stageReady;
            if (token !== this.renderToken) {
                return;
            }
            const queryAlignedColor = toMolstarColor(this.queryAlignedColor, DEFAULT_QUERY_COLOR);
            const queryUnalignedColor = toMolstarColor(this.queryUnalignedColor, DEFAULT_QUERY_UNALIGNED_COLOR);
            await this.renderQueryVisuals(token, queryAlignedColor, queryUnalignedColor, this.isMultimerMode());
        },
        buildHighlightExpression() {
            if (!Array.isArray(this.highlights) || this.highlights.length === 0) {
                return null;
            }
            const chainSet = this.targetChains;
            const expressions = [];
            this.highlights.forEach((highlight, index) => {
                if (!highlight || highlight.length < 2) return;
                const start = Number(highlight[0]);
                const length = Number(highlight[1]);
                if (!Number.isFinite(start) || !Number.isFinite(length) || length <= 0) return;
                const alignment = this.alignments[index];
                if (!alignment) return;
                const chain = this.getAlignmentChain('target', alignment, index);
                let chainMap = this.targetSerialMap.get(chain);
                if (!chainMap && this.targetSerialMap.size === 1) {
                    chainMap = Array.from(this.targetSerialMap.values())[0];
                }
                const mappedRanges = mapRangesToAuth([{ start, end: start + length - 1 }], chainMap);
                if (!mappedRanges || mappedRanges.length === 0) return;
                const useChain = chainSet && chainSet.size > 0 && chainSet.has(chain);
                expressions.push(buildChainExpression(chain, mappedRanges, useChain));
            });
            if (expressions.length === 0) return null;
            return expressions.length === 1
                ? expressions[0]
                : MS.struct.combinator.merge(expressions);
        },
        focusSelection() {
            if (!this.stage) return;
            const selectionTarget = this.queryStructureRef || this.targetStructureRef;
            if (!selectionTarget) {
                this.stage.focusLoci(null, this.autoViewTime);
                return;
            }
            const expr = this.queryStructureRef
                ? this.buildSelectionExpression('query', this.showQuery)
                : this.buildSelectionExpression('target', this.showTarget);
            const loci = getSelectionLoci(expr, selectionTarget);
            this.stage.focusLoci(loci, this.autoViewTime);
        },
        async loadStructureData() {
            if (!this.alignments || this.alignments.length === 0) {
                return;
            }

            this.updateAlignmentChains();
            this.queryData = null;
            this.targetData = null;
            this.lastQueryPdb = '';
            this.lastTargetPdb = '';
            this.hasQuery = true;
            this.tmAlignResults = null;

            let queryPdb = "";
            if (this.$LOCAL) {
                if (this.hits.queries[0].hasOwnProperty('pdb')) {
                    queryPdb = JSON.parse(this.hits.queries[0].pdb);
                } else {
                    queryPdb = mockPDB(this.hits.queries[0].qCa, this.hits.queries[0].sequence, 'A');
                }
            } else if (this.$route.params.ticket.startsWith('user-')) {
                if (this.hits.queries[0].hasOwnProperty('pdb')) {
                    queryPdb = JSON.parse(this.hits.queries[0].pdb);
                } else {
                    const localData = this.$root.userData[this.$route.params.entry];
                    queryPdb = mockPDB(localData.queries[0].qCa, localData.queries[0].sequence, 'A');
                }
            } else {
                try {
                    const request = await this.$axios.get("api/result/" + this.$route.params.ticket + '/query');
                    queryPdb = request.data;
                } catch (e) {
                    queryPdb = "";
                    this.hasQuery = false;
                }
            }

            if (queryPdb) {
                const trimmed = queryPdb.trimStart();
                const format = (trimmed.startsWith('#') || trimmed.startsWith('data_')) ? 'mmcif' : 'pdb';
                this.queryData = { data: queryPdb, format, label: 'query' };
                this.lastQueryPdb = queryPdb;
                const parsed = await readQueryDataWithMolIo(queryPdb, format);
                if (parsed) {
                    this.queryChains = parsed.chainIds;
                    this.querySerialMap = parsed.serialMap;
                    this.updateQueryAlignmentChainsFromStructure();
                } else {
                    this.queryChains = new Set();
                    this.querySerialMap = new Map();
                }
            } else {
                this.hasQuery = false;
                this.queryChains = new Set();
                this.querySerialMap = new Map();
            }

            const targets = [];
            let lastIdx = null;
            let remoteData = null;
            let i = 0;
            for (let idx = 0; idx < this.alignments.length; idx += 1) {
                const alignment = this.alignments[idx];
                const chain = this.getAlignmentChain('target', alignment, idx);
                let tSeq = alignment.tSeq;
                let tCa = alignment.tCa;
                if (Number.isInteger(alignment.tCa) && Number.isInteger(alignment.tSeq)) {
                    const db = alignment.db;
                    const remoteIdx = alignment.tCa;
                    if (remoteIdx != lastIdx) {
                        const ticket =  this.$route.params.ticket;
                        const response = await this.$axios.get("api/result/" + ticket + '/' + this.$route.params.entry + '?format=brief&index=' + remoteIdx + '&database=' + db);
                        remoteData = response.data;
                        lastIdx = remoteIdx;
                    }
                    // console.log("remoteData", remoteData, lastIdx)
                    if (Array.isArray(remoteData) && remoteData.length > 0) {
                        tSeq = remoteData[i].tSeq;
                        tCa = remoteData[i].tCa;
                        i++;
                    } else {
                        tSeq = null;
                        tCa = null;
                    }
                }
                // console.log(tSeq, tCa)
                if (!tSeq || !tCa) {
                    continue;
                }
                const mock = mockPDB(tCa, tSeq, chain);
                try {
                    const pdb = await pulchra(mock);
                    targets.push(applyChainId(pdb, chain));
                } catch (error) {
                    console.warn('pulchra failed for target chain', chain, error);
                    targets.push(applyChainId(mock, chain));
                }
            }
            let targetPdb = mergePdbChunks(targets);
            const isMultimer = this.isMultimerMode();
            if (isMultimer) {
                const complexT = typeof this.alignments[0].complext === 'string' ? this.alignments[0].complext : '';
                const complexU = typeof this.alignments[0].complexu === 'string' ? this.alignments[0].complexu : '';
                const t = complexT.split(',').map(x => parseFloat(x));
                const uFlat = complexU.split(',').map(x => parseFloat(x));
                if (t.length === 3 && uFlat.length === 9 && t.every(Number.isFinite) && uFlat.every(Number.isFinite)) {
                    const u = [
                        [uFlat[0], uFlat[1], uFlat[2]],
                        [uFlat[3], uFlat[4], uFlat[5]],
                        [uFlat[6], uFlat[7], uFlat[8]],
                    ];
                    targetPdb = transformPdb(targetPdb, t, u);
                }
            }
            if (targetPdb) {
                this.targetChains = getChainIdsFromPdb(targetPdb);
                this.targetSerialMap = buildSerialResidueMap(targetPdb);
                this.updateTargetAlignmentChainsFromStructure();
                if (this.hasQuery && this.queryData && !isMultimer) {
                    if (this.queryData.format === 'pdb') {
                        try {
                            const alnFasta = `>target\n${this.alignments[0].dbAln}\n\n>query\n${this.alignments[0].qAln}`;
                            const tRanges = buildAlignmentRanges(this.alignments, 'target', (alignment, index) => this.getAlignmentChain('target', alignment, index));
                            const qRanges = buildAlignmentRanges(this.alignments, 'query', (alignment, index) => this.getAlignmentChain('query', alignment, index));
                            const tSubPdb = makeSubPdbFromRanges(targetPdb, tRanges);
                            const qSubPdb = makeSubPdbFromRanges(this.lastQueryPdb, qRanges);
                            // console.log(targetPdb)
                            // console.log(this.lastQueryPdb)
                            // console.log("target", tRanges, targetPdb)
                            // console.log("query", qRanges, this.lastQueryPdb)
                            const tm = await tmalign(tSubPdb, qSubPdb, alnFasta);
                            this.tmAlignResults = parseTMOutput(tm.output);
                            const { t, u } = parseTMMatrix(tm.matrix);
                            targetPdb = transformPdb(targetPdb, t, u);
                        } catch (error) {
                            console.warn('tmalign-wasm failed, using unaligned target', error);
                        }
                    } else {
                        console.warn('Skipping TM-align: query structure is not PDB');
                    }
                }
                this.targetData = { data: targetPdb, format: 'pdb', label: 'target' };
                this.lastTargetPdb = targetPdb;
            } else {
                this.targetChains = new Set();
                this.targetSerialMap = new Map();
                this.targetAlignmentChains = new Map();
            }
        },
        async renderStructures(focus = true, reloadStructures = true) {
            if (!this.stage || !this.stageReady) {
                return;
            }
            const token = ++this.renderToken;
            await this.stageReady;
            if (token !== this.renderToken) {
                return;
            }

            const queryAlignedColor = toMolstarColor(this.queryAlignedColor, DEFAULT_QUERY_COLOR);
            const queryUnalignedColor = toMolstarColor(this.queryUnalignedColor, DEFAULT_QUERY_UNALIGNED_COLOR);
            const targetAlignedColor = toMolstarColor(this.targetAlignedColor, DEFAULT_TARGET_COLOR);
            const targetUnalignedColor = toMolstarColor(this.targetUnalignedColor, DEFAULT_TARGET_UNALIGNED_COLOR);
            const highlightColor = toMolstarColor(DEFAULT_HIGHLIGHT_COLOR, DEFAULT_HIGHLIGHT_COLOR);
            const isMultimer = this.isMultimerMode();

            const needLoadTarget = reloadStructures || (this.targetData && !this.targetStructureRef);
            const needLoadQuery = reloadStructures || (this.queryData && !this.queryStructureRef);
            if (reloadStructures || needLoadTarget || needLoadQuery) {
                await this.stage.clear();
                this.queryStructureRef = null;
                this.targetStructureRef = null;
                this.queryVisualComponentRefs = [];
                this.targetVisualComponentRefs = [];
                this.highlightComponentRef = null;
                this.arrowRefs = [];
                if (this.targetData) {
                    this.targetStructureRef = await this.stage.loadStructure(this.targetData);
                }
                if (this.queryData) {
                    this.queryStructureRef = await this.stage.loadStructure(this.queryData);
                }
            } else {
                await this.clearVisualComponents();
            }

            await this.renderTargetVisuals(token, targetAlignedColor, targetUnalignedColor, isMultimer);
            await this.renderQueryVisuals(token, queryAlignedColor, queryUnalignedColor, isMultimer);

            await this.rerenderHighlight(highlightColor);
            await this.rerenderArrows();

            const focusRef = this.queryStructureRef || this.targetStructureRef;
            if (focus) {
                const loci = getSelectionLoci(focusRef ? MS.struct.generator.all() : null, focusRef);
                this.stage.focusLoci(loci, this.autoViewTime);
            }
        },
    },
    watch: {
        'showQuery': function() {
            if (!this.stage) {
                return;
            }
            this.rerenderQueryVisuals();
        },
        'showTarget': function(val, _) {
            if (!this.stage) {
                return;
            }
            this.rerenderTargetVisuals();
        },
        'showArrows': function() {
            if (!this.stage) {
                return;
            }
            this.rerenderArrows();
        },
        'highlights': {
            deep: true,
            handler() {
                if (!this.stage) {
                    return;
                }
                this.rerenderHighlight();
            },
        }
    },
    computed: {
        querySele: function() {
            if (this.alignments.length === 0 || this.showQuery == 2)
                return '';
            if (this.showQuery === 0)
                return this.alignments.map((a, index) => `${a.qStartPos}-${a.qEndPos}:${this.getAlignmentChain('query', a, index)}`).join(" or ");
            if (this.showQuery === 1)
                return this.alignments.map((a, index) => `:${this.getAlignmentChain('query', a, index)}`).join(" or ");
        },
        targetSele: function() {
            if (this.alignments.length === 0 || this.showTarget == 2)
                return '';
            if (this.showTarget === 0)
                return this.alignments.map((a, index) => `${a.dbStartPos}-${a.dbEndPos}:${this.getAlignmentChain('target', a, index)}`).join(" or ");
            if (this.showTarget === 1)
                return this.alignments.map((a, index) => `:${this.getAlignmentChain('target', a, index)}`).join(" or ");
        },
        tmPanelBindings: function() {
            return (this.isFullscreen) ? { 'style': 'margin-top: 10px; font-size: 2em; line-height: 2em' } : { }
        },
    },
    async mounted() {
        if (!this.alignments || this.alignments.length === 0) {
            return;
        }
        if (typeof(this.alignments[0].tCa) == "undefined") {
            return;
        }
        await this.stageReady;
        await this.loadStructureData();
        await this.renderStructures(true, true);
    },
}
</script>

<style scoped>
.structure-wrapper {
    width: 500px;
    height: 400px;
    margin: 0 auto;
}
</style>

<style>
.theme--dark .structure-wrapper .v-tooltip__content {
    background: rgba(97, 97, 97, 0.3);
}
/* @media only screen and (max-width: 600px) {
    .structure-wrapper {
        width: 300px;
    }
} */
.structure-viewer {
    width: 100%;
    height: 100%;
}
.structure-viewer canvas {
    border-radius: 2px;
}
.structure-panel {
    position: relative;
}
.toolbar-panel {
    display: inline-flex;
    flex-direction: row;
    position: absolute;
    justify-content: center;
    width: 100%;
    bottom: 0;
    z-index: 1;
    left: 0;
}
.tmscore-panel {
    position: absolute;
    width: 100%;
    top: 0;
    left: 0;
    z-index: 1;
    font-family: monospace;
    color: rgb(31, 119, 180);
}
.left-cell {
    text-align: right;
    width: 50%;
}
.right-cell {
    text-align: left;
    width: 50%;
    padding-left: 0.3em;
}
</style>
