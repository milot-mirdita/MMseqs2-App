<template>
<div class="structure-panel">
    <StructureViewerTooltip attach=".structure-panel" />
    <div class="structure-wrapper" ref="structurepanel">
        <StructureViewerToolbar
            :isFullscreen="isFullscreen"
            :isSpinning="isSpinning"
            @makeImage="handleMakeImage"
            @makePDB="handleMakePDB"
            @resetView="handleResetView"
            @toggleFullscreen="handleToggleFullscreen"
            @toggleSpin="handleToggleSpin"
            disableArrowButton
            disableQueryButton
            disableTargetButton
            style="position: absolute; bottom: 8px;"
        />
        <div class="structure-viewer" ref="viewport"></div>
    </div>
</div>
</template>

<script>
import StructureViewerTooltip from './StructureViewerTooltip.vue';
import StructureViewerToolbar from './StructureViewerToolbar.vue';
import StructureViewerMixin from './StructureViewerMixin.vue';
import { tmalign, parse as parseTMOutput, parseMatrix as parseTMMatrix } from 'tmalign-wasm';
import { mockPDB, downloadBlob } from './Utilities.js';
import { pulchra } from 'pulchra-wasm';
import { MolScriptBuilder as MS } from 'molstar/lib/mol-script/language/builder';
import { Script } from 'molstar/lib/mol-script/script';
import { StructureSelection, StructureElement, StructureProperties } from 'molstar/lib/mol-model/structure';
import { Color } from 'molstar/lib/mol-util/color';

const DEFAULT_REFERENCE_COLOR = 0x1e88e5;
const DEFAULT_REGULAR_COLOR = 0xffc107;
const DEFAULT_MASK_COLOR = 0x666666;
const DEFAULT_HIGHLIGHT_COLOR = 0x11ffee;

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
    if (!pdb) return map;
    const lastResidue = new Map();
    for (const line of pdb.split('\n')) {
        const parsed = parseAtomLine(line);
        if (!parsed) continue;
        const chainKey = parsed.chainId && parsed.chainId.trim() ? parsed.chainId.trim() : '_';
        const residueKey = `${chainKey}:${parsed.resno}:${parsed.resname}`;
        if (lastResidue.get(chainKey) === residueKey) continue;
        lastResidue.set(chainKey, residueKey);
        if (!map.has(chainKey)) map.set(chainKey, []);
        map.get(chainKey).push(parsed.resno);
    }
    return map;
};

const buildSerialIndexMap = (serialMap) => {
    const map = new Map();
    if (!serialMap) return map;
    serialMap.forEach((resno, idx) => {
        if (!map.has(resno)) {
            map.set(resno, idx);
        }
    });
    return map;
};

const mapRangesToAuth = (ranges, serialMap) => {
    if (!serialMap || serialMap.length === 0) return ranges;
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

const positionsToRanges = (positions) => {
    if (!positions || positions.length === 0) return [];
    const sorted = Array.from(new Set(positions)).sort((a, b) => a - b);
    const ranges = [];
    let start = sorted[0];
    let end = start;
    for (let i = 1; i < sorted.length; i += 1) {
        const pos = sorted[i];
        if (pos === end + 1) {
            end = pos;
            continue;
        }
        ranges.push({ start: start + 1, end: end + 1 });
        start = pos;
        end = pos;
    }
    ranges.push({ start: start + 1, end: end + 1 });
    return ranges;
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
    if (!data || !expression) return null;
    const selection = Script.getStructureSelection(expression, data);
    return StructureSelection.toLociWithSourceUnits(selection);
};

const extractAtomLines = (pdb) => {
    if (!pdb) return [];
    return pdb.split('\n').filter(line => line.startsWith('ATOM'));
};

const transformPdb = (pdb, t, u) => {
    if (!pdb) return pdb;
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

const getChainIdsFromPdb = (pdb) => {
    const chains = new Set();
    if (!pdb) return chains;
    for (const line of pdb.split('\n')) {
        const parsed = parseAtomLine(line);
        if (!parsed) continue;
        if (parsed.chainId && parsed.chainId.trim()) {
            chains.add(parsed.chainId.trim());
        }
    }
    return chains;
};

const makeSubPdbFromRanges = (pdb, rangesByChain) => {
    if (!pdb) return '';
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

// Mock alignment object from two (MSA-derived) aligned strings
function mockAlignment(one, two) {
    let res = { backtrace: "", qAln: "", dbAln: "" };
    let started = false; // flag for first Match column in backtrace
    let m = 0;           // index in msa
    let qr = 0;          // index in seq
    let tr = 0;
    let qBuffer = "";
    let tBuffer = "";
    while (m < one.length) {
        const qc = one[m];
        const tc = two[m];
        if (qc === '-' && tc === '-') {
            // Skip gap columns
        } else if (qc === '-') {
            if (started) {
                res.backtrace += 'D';               
                qBuffer += qc;
                tBuffer += tc;
            }
            ++tr;
        } else if (tc === '-') {
            if (started) {
                res.backtrace += 'I';
                qBuffer += qc;
                tBuffer += tc;
            }
            ++qr;
        } else {
            if (started) {
                res.qAln += qBuffer;
                res.dbAln += tBuffer;
                qBuffer = "";
                tBuffer = "";
            } else {
                started = true;
                res.qStartPos = qr;
                res.dbStartPos = tr;
            }
            res.backtrace += 'M';
            qBuffer += qc;
            tBuffer += tc;
            res.qEndPos = qr;
            res.dbEndPos = tr;
            ++qr;
            ++tr;
        }
        ++m;
    }
    res.qStartPos++;
    res.dbStartPos++;
    res.qSeq  = one.replace(/-/g, '');
    res.tSeq  = two.replace(/-/g, '');
    return res;
}


function getMaskedPositions(seq, mask) {
    const result = [];
    let resno = 0;
    for (let i = 0; i < seq.length; i++) {
        if (seq[i] !== '-') {
            if (mask[i] === 0) {
                result.push(resno);
            }
            resno++;
        }
    }
    return result;
}

function getAlignmentPos(seq, residueIndex) {
    let resno = -1;
    for (let i = 0; i < seq.length; i++) {
        if (seq[i] !== '-') {
            resno++;
        }
        if (resno == residueIndex) {
            return i;
        }
    }
    return -1;
}

function getResidueIndex(seq, alignmentPos) {
    let residueIndex = -1;
    for (let i = 0; i <= alignmentPos && i < seq.length; i++) {
        if (seq[i] !== '-') {
            residueIndex++;
        }
    }
    return residueIndex;
}

export default {
    name: "StructureViewerMSA",
    components: {
        StructureViewerToolbar,
        StructureViewerTooltip,
    },
    mixins: [
        StructureViewerMixin,
    ],
    data: () => ({
        structureItems: [],
        structureIndexByStructure: new Map(),
        renderToken: 0,
        overlayToken: 0,
        clickUnsub: null,
        pdbCache: new Map(),
        selectedColumn: -1,
    }),
    props: {
        entries: { type: Array, required: true },
        selection: { type: Array, required: true, default: [0, 1] },
        mask: { type: Array, required: true },
        reference: { type: Number, required: true },
        bgColorLight: { type: String, default: "white" },
        bgColorDark: { type: String, default: "#1E1E1E" },
        representationStyle: { type: String, default: "cartoon" },
        referenceStyleParameters: {
            type: Object,
            default: () => ({ color: 0x1E88E5, opacity: 1.0 })
        },
        regularStyleParameters: {
            type: Object,
            default: () => ({ color: 0xFFC107, opacity: 0.5, side: 'front' })
        },
    },
    async mounted() {
        await this.stageReady;
        this.subscribeClicks();
        await this.rebuildStructures(true);
    },
    beforeDestroy() {
        if (this.clickUnsub) {
            this.clickUnsub();
            this.clickUnsub = null;
        }
    },
    methods: {
        resetView() {
            if (!this.stage) return;
            this.focusReference();
        },
        makePDB() {
            if (!this.structureItems.length) return;
            let result = `\
TITLE     Superposed structures from Foldmason alignment
REMARK    This file was generated by the FoldMason webserver:
REMARK      https://search.foldseek.com/foldmason
REMARK    Please cite:
REMARK      https://doi.org/10.1101/2024.08.01.606130
REMARK    Warning: Non C-alpha atoms may have been re-generated by PULCHRA
REMARK             if they are not present in the original PDB file.
`;
            for (const item of this.structureItems) {
                const entry = this.entries[item.index];
                const PDB = extractAtomLines(item.pdb).join('\n');
                const name = entry?.name || `key-${item.index}`;
                const remark = `REMARK    Name: ${name}`;
                result += `\
MODEL     ${item.index}
${remark}
${PDB}
ENDMDL
`;
            }
            result += "END";
            downloadBlob(new Blob([result], { type: 'text/plain' }), "foldmason.pdb");
        },
        async makeImage() {
            if (!this.stage) return;
            const wasSpinning = this.isSpinning;
            this.isSpinning = false;
            const blob = await this.stage.makeImage();
            if (blob) {
                downloadBlob(blob, "foldmason.png");
            }
            this.isSpinning = wasSpinning;
        },
        subscribeClicks() {
            if (this.clickUnsub) {
                this.clickUnsub();
            }
            this.clickUnsub = this.stage.onClick((event) => this.handleStructureClick(event));
        },
        handleStructureClick(event) {
            const loci = event?.current?.loci || event?.current;
            if (!StructureElement.Loci.is(loci)) {
                this.selectedColumn = -1;
                this.$emit('columnSelected', -1);
                this.renderOverlays();
                return;
            }
            const location = StructureElement.Loci.getFirstLocation(loci);
            if (!location) {
                this.selectedColumn = -1;
                this.$emit('columnSelected', -1);
                this.renderOverlays();
                return;
            }
            const structure = location.structure;
            const index = this.structureIndexByStructure.get(structure);
            if (index === undefined) {
                return;
            }
            const chainId = StructureProperties.chain.auth_asym_id(location);
            const resno = StructureProperties.residue.auth_seq_id(location);
            const item = this.structureItems.find(entry => entry.index === index);
            const serialIndex = item ? this.getSerialIndex(item, chainId, resno) : null;
            if (!Number.isFinite(serialIndex)) {
                this.selectedColumn = -1;
                this.$emit('columnSelected', -1);
                this.renderOverlays();
                return;
            }
            const alnPos = getAlignmentPos(this.entries[index].aa, serialIndex);
            this.selectedColumn = alnPos;
            this.$emit('columnSelected', alnPos);
            this.renderOverlays();
        },
        focusReference() {
            if (!this.stage) return;
            const refItem = this.structureItems.find(item => item.index === this.reference) || this.structureItems[0];
            if (!refItem) return;
            const loci = getSelectionLoci(MS.struct.generator.all(), refItem.structureRef);
            this.stage.focusLoci(loci, this.transitionDuration);
        },
        getPrimaryChain(item) {
            if (item.primaryChain) return item.primaryChain;
            if (item.chainIds && item.chainIds.size > 0) {
                return Array.from(item.chainIds.values())[0];
            }
            return 'A';
        },
        getChainMap(item, chain) {
            if (!item.serialMap || item.serialMap.size === 0) return null;
            let chainMap = item.serialMap.get(chain);
            if (!chainMap && item.serialMap.size === 1) {
                chainMap = Array.from(item.serialMap.values())[0];
            }
            return chainMap || null;
        },
        getSerialIndex(item, chain, resno) {
            if (!item.serialIndexByChain || item.serialIndexByChain.size === 0) return null;
            let chainMap = item.serialIndexByChain.get(chain);
            if (!chainMap && item.serialIndexByChain.size === 1) {
                chainMap = Array.from(item.serialIndexByChain.values())[0];
            }
            if (!chainMap) return null;
            const value = chainMap.get(resno);
            return Number.isFinite(value) ? value : null;
        },
        buildMaskExpression(entry, item) {
            if (!this.mask || this.mask.length === 0) return null;
            const positions = getMaskedPositions(entry.aa, this.mask);
            if (!positions.length) return null;
            const ranges = positionsToRanges(positions);
            const chain = this.getPrimaryChain(item);
            const chainMap = this.getChainMap(item, chain);
            const mappedRanges = mapRangesToAuth(ranges, chainMap);
            if (!mappedRanges.length) return null;
            const useChain = item.chainIds && item.chainIds.size > 0 && item.chainIds.has(chain);
            return buildChainExpression(chain, mappedRanges, useChain);
        },
        buildHighlightExpression(entry, item) {
            if (this.selectedColumn < 0) return null;
            const residueIndex = getResidueIndex(entry.aa, this.selectedColumn);
            if (!Number.isFinite(residueIndex) || residueIndex < 0) return null;
            const chain = this.getPrimaryChain(item);
            const chainMap = this.getChainMap(item, chain);
            const mappedRanges = mapRangesToAuth([{ start: residueIndex + 1, end: residueIndex + 1 }], chainMap);
            if (!mappedRanges.length) return null;
            const useChain = item.chainIds && item.chainIds.size > 0 && item.chainIds.has(chain);
            return buildChainExpression(chain, mappedRanges, useChain);
        },
        async buildEntryPdb(index, entry) {
            if (this.pdbCache.has(index)) {
                return this.pdbCache.get(index);
            }
            const seq = entry.aa ? entry.aa.replace(/-/g, '') : '';
            const mock = mockPDB(entry.ca, seq, 'A');
            let pdb = mock;
            try {
                pdb = await pulchra(mock);
            } catch (error) {
                console.warn('pulchra failed for entry', entry?.name || index, error);
            }
            this.pdbCache.set(index, pdb);
            return pdb;
        },
        async alignEntryToReference(refEntry, refPdb, entry, entryPdb) {
            const aln = mockAlignment(refEntry.aa, entry.aa);
            if (!Number.isFinite(aln.qStartPos) || !Number.isFinite(aln.dbStartPos)) {
                return entryPdb;
            }
            const refChains = getChainIdsFromPdb(refPdb);
            const entryChains = getChainIdsFromPdb(entryPdb);
            const refChain = refChains.size ? Array.from(refChains)[0] : 'A';
            const entryChain = entryChains.size ? Array.from(entryChains)[0] : 'A';
            const qRanges = new Map([[refChain, [{ start: aln.qStartPos, end: aln.qEndPos }]]]);
            const tRanges = new Map([[entryChain, [{ start: aln.dbStartPos, end: aln.dbEndPos }]]]);
            const qSubPdb = makeSubPdbFromRanges(refPdb, qRanges);
            const tSubPdb = makeSubPdbFromRanges(entryPdb, tRanges);
            const alnFasta = `>target\n${aln.dbAln}\n\n>query\n${aln.qAln}`;
            try {
                const tm = await tmalign(tSubPdb, qSubPdb, alnFasta);
                const { t, u } = parseTMMatrix(tm.matrix);
                parseTMOutput(tm.output);
                return transformPdb(entryPdb, t, u);
            } catch (error) {
                console.warn('tmalign-wasm failed for entry', entry?.name || '', error);
            }
            return entryPdb;
        },
        async addStructureItem(index, pdb, isReference) {
            if (!this.stage) return null;
            const structureRef = await this.stage.loadStructure({ data: pdb, format: 'pdb', label: `key-${index}` });
            if (!structureRef) return null;
            const serialMap = buildSerialResidueMap(pdb);
            const serialIndexByChain = new Map();
            serialMap.forEach((values, chain) => {
                serialIndexByChain.set(chain, buildSerialIndexMap(values));
            });
            const chainIds = getChainIdsFromPdb(pdb);
            const primaryChain = chainIds.size ? Array.from(chainIds.values())[0] : 'A';
            const baseComponent = await this.stage.createComponentStatic(structureRef, 'all');
            const baseColor = toMolstarColor(
                isReference ? this.referenceStyleParameters.color : this.regularStyleParameters.color,
                isReference ? DEFAULT_REFERENCE_COLOR : DEFAULT_REGULAR_COLOR
            );
            const opacityValue = isReference ? this.referenceStyleParameters.opacity : this.regularStyleParameters.opacity;
            const initialState = Number.isFinite(opacityValue) ? { alphaFactor: opacityValue } : undefined;
            if (baseComponent) {
                await this.stage.addRepresentation(
                    baseComponent,
                    {
                        type: this.representationStyle,
                        color: 'uniform',
                        colorParams: { value: baseColor },
                    },
                    initialState ? { initialState } : undefined
                );
            }
            const item = {
                index,
                pdb,
                structureRef,
                chainIds,
                primaryChain,
                serialMap,
                serialIndexByChain,
                baseComponent,
                overlays: {},
            };
            this.structureItems.push(item);
            const structure = structureRef?.cell?.obj?.data;
            if (structure) {
                this.structureIndexByStructure.set(structure, index);
            }
            return item;
        },
        async renderOverlays() {
            if (!this.stage || !this.stageReady) return;
            const token = ++this.overlayToken;
            await this.stageReady;
            if (token !== this.overlayToken) return;
            for (const item of this.structureItems) {
                if (item.overlays?.mask) {
                    await this.stage.remove(item.overlays.mask);
                }
                if (item.overlays?.highlight) {
                    await this.stage.remove(item.overlays.highlight);
                }
                item.overlays = {};
            }
            if (token !== this.overlayToken) return;
            for (const item of this.structureItems) {
                const entry = this.entries[item.index];
                if (!entry) continue;
                const maskExpr = this.buildMaskExpression(entry, item);
                if (maskExpr) {
                    const maskComponent = await this.stage.createComponentFromExpression(
                        item.structureRef,
                        maskExpr,
                        `mask-${item.index}-${token}`
                    );
                    if (maskComponent) {
                        await this.stage.addRepresentation(maskComponent, {
                            type: this.representationStyle,
                            color: 'uniform',
                            colorParams: { value: Color(DEFAULT_MASK_COLOR) },
                        });
                        item.overlays.mask = maskComponent;
                    }
                }
                const highlightExpr = this.buildHighlightExpression(entry, item);
                if (highlightExpr) {
                    const highlightComponent = await this.stage.createComponentFromExpression(
                        item.structureRef,
                        highlightExpr,
                        `highlight-${item.index}-${token}`
                    );
                    if (highlightComponent) {
                        await this.stage.addRepresentation(highlightComponent, {
                            type: 'ball-and-stick',
                            color: 'uniform',
                            colorParams: { value: Color(DEFAULT_HIGHLIGHT_COLOR) },
                        });
                        item.overlays.highlight = highlightComponent;
                    }
                }
            }
        },
        async rebuildStructures(focus = true) {
            if (!this.stage || !this.stageReady) return;
            const token = ++this.renderToken;
            await this.stageReady;
            if (token !== this.renderToken) return;

            await this.stage.clear();
            this.structureItems = [];
            this.structureIndexByStructure = new Map();

            const indices = [];
            if (Number.isInteger(this.reference)) {
                indices.push(this.reference);
            }
            if (Array.isArray(this.selection)) {
                for (const idx of this.selection) {
                    if (idx !== this.reference) {
                        indices.push(idx);
                    }
                }
            }

            if (indices.length === 0) return;
            const refEntry = this.entries[this.reference];
            if (!refEntry) return;

            const refPdb = await this.buildEntryPdb(this.reference, refEntry);
            await this.addStructureItem(this.reference, refPdb, true);

            for (const idx of indices) {
                if (idx === this.reference) continue;
                const entry = this.entries[idx];
                if (!entry) continue;
                const entryPdb = await this.buildEntryPdb(idx, entry);
                const alignedPdb = await this.alignEntryToReference(refEntry, refPdb, entry, entryPdb);
                await this.addStructureItem(idx, alignedPdb, false);
            }

            await this.renderOverlays();

            if (focus) {
                this.focusReference();
            }
        },
    },
    watch: {
        selection: function(newV, oldV) {
            this.rebuildStructures(true);
        },
        mask: function(newM, oldM) {
            this.renderOverlays();
        },
        reference: function() {
            this.rebuildStructures(true);
        },
        entries: {
            deep: true,
            handler() {
                this.pdbCache.clear();
                this.rebuildStructures(true);
            },
        },
    },
}
</script>

<style scoped>
.structure-panel {
    width: 100%;
    height: 100%;
    position: relative;
}
.structure-viewer {
    width: 100%;
    height: 100%;
    position: absolute;
    bottom: 0;
    left: 0;
    right: 0;
    top: 0;
    margin: 0;
    padding: 0;
    overflow: hidden;
}
</style>
