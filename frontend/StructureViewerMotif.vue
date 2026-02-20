<template>
    <div class="structure-panel">
        <StructureViewerTooltip attach=".structure-panel" />
        <div class="structure-wrapper" ref="structurepanel">
            <table class="tmscore-panel" v-bind="tmPanelBindings">
                <tr>
                    <td class="left-cell">idf-score:</td>
                    <td class="right-cell">{{ alignment.idfscore }}</td>
                </tr>
                <tr>
                    <td class="left-cell">RMSD:</td>
                    <td class="right-cell">{{ alignment.rmsd }}</td>
                </tr>
            </table>
            <StructureViewerToolbar
                :isFullscreen="isFullscreen"
                :isSpinning="isSpinning"
                :showTarget="showTarget"
                :showQuery="showQuery"
                :disableArrowButton="true"
                @makeImage="handleMakeImage"
                @makePDB="handleMakePDB"
                @resetView="handleResetView"
                @toggleFullscreen="handleToggleFullscreen"
                @toggleQuery="handleToggleQuery"
                @toggleTarget="handleToggleTarget"
                @toggleSpin="handleToggleSpin"
                style="position: absolute; bottom: 8px;"
            />
            <div class="structure-viewer" ref="viewport"></div>
        </div>
    </div>
</template>

<script>
import StructureViewerMixin from './StructureViewerMixin.vue';
import StructureViewerToolbar from './StructureViewerToolbar.vue';
import StructureViewerTooltip from './StructureViewerTooltip.vue';
import { downloadBlob, splitAlphaNum } from './Utilities.js';
import { MolScriptBuilder as MS } from 'molstar/lib/mol-script/language/builder';
import { Script } from 'molstar/lib/mol-script/script';
import { StructureSelection } from 'molstar/lib/mol-model/structure';
import { Color } from 'molstar/lib/mol-util/color';

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

const normalizeRepresentation = (repr, fallback) => {
    const value = (repr || '').toLowerCase();
    if (value === 'ribbon' || value === 'cartoon' || value === 'tube') return 'cartoon';
    if (value === 'licorice' || value === 'ball-and-stick') return 'ball-and-stick';
    if (value === 'surface') return 'molecular-surface';
    return fallback;
};

const processPdb = (rawpdb) => {
    if (!rawpdb || typeof rawpdb !== 'string') {
        return null;
    }
    const trimmed = rawpdb.trimStart();
    if (!trimmed) {
        return null;
    }
    if (trimmed[0] === '#' || trimmed.startsWith('data_')) {
        return { data: trimmed, format: 'mmcif' };
    }
    return { data: rawpdb, format: 'pdb' };
};

const extractAtomLines = (pdb) => {
    if (!pdb) return [];
    return pdb.split('\n').filter(line => line.startsWith('ATOM') || line.startsWith('HETATM'));
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

const motifToMap = (motif) => {
    const byChain = new Map();
    if (!motif || typeof motif !== 'string') {
        return byChain;
    }
    motif.split(',').forEach(token => {
        const part = token.trim();
        if (!part) return;
        const core = part.split(':')[0];
        const [chain, pos] = splitAlphaNum(core);
        const resno = Number.parseInt(pos, 10);
        if (!chain || !Number.isFinite(resno)) return;
        if (!byChain.has(chain)) byChain.set(chain, new Set());
        byChain.get(chain).add(resno);
    });
    return byChain;
};

const residuesToRanges = (residueSet) => {
    const residues = Array.from(residueSet).sort((a, b) => a - b);
    if (residues.length === 0) return [];
    const ranges = [];
    let start = residues[0];
    let end = start;
    for (let i = 1; i < residues.length; i += 1) {
        const value = residues[i];
        if (value === end + 1) {
            end = value;
            continue;
        }
        ranges.push({ start, end });
        start = value;
        end = value;
    }
    ranges.push({ start, end });
    return ranges;
};

const buildChainExpression = (chain, ranges = null) => {
    const groupBy = MS.struct.atomProperty.macromolecular.residueKey();
    const chainTest = MS.core.rel.eq([
        MS.struct.atomProperty.macromolecular.auth_asym_id(),
        chain,
    ]);
    if (!ranges || ranges.length === 0) {
        return MS.struct.generator.atomGroups({ 'chain-test': chainTest, 'group-by': groupBy });
    }
    const parts = ranges.map(range => MS.struct.generator.atomGroups({
        'chain-test': chainTest,
        'residue-test': MS.core.rel.inRange([
            MS.struct.atomProperty.macromolecular.auth_seq_id(),
            range.start,
            range.end,
        ]),
        'group-by': groupBy,
    }));
    return parts.length === 1 ? parts[0] : MS.struct.combinator.merge(parts);
};

const buildMotifExpression = (motif) => {
    const byChain = motifToMap(motif);
    const expressions = [];
    byChain.forEach((residues, chain) => {
        const ranges = residuesToRanges(residues);
        expressions.push(buildChainExpression(chain, ranges));
    });
    if (expressions.length === 0) return null;
    return expressions.length === 1 ? expressions[0] : MS.struct.combinator.merge(expressions);
};

const buildMotifChainExpression = (motif) => {
    const byChain = motifToMap(motif);
    const expressions = [];
    byChain.forEach((_, chain) => {
        expressions.push(buildChainExpression(chain));
    });
    if (expressions.length === 0) return null;
    return expressions.length === 1 ? expressions[0] : MS.struct.combinator.merge(expressions);
};

const getSelectionLoci = (expression, structureRef) => {
    const data = structureRef?.cell?.obj?.data;
    if (!data || !expression) return null;
    const selection = Script.getStructureSelection(expression, data);
    return StructureSelection.toLociWithSourceUnits(selection);
};

export default {
    name: 'StructureViewerMotif',
    components: {
        StructureViewerToolbar,
        StructureViewerTooltip,
    },
    mixins: [StructureViewerMixin],
    data: () => ({
        showQuery: 0,
        showTarget: 0,
        renderToken: 0,
        queryStructureRef: null,
        targetStructureRef: null,
        lastQueryPdb: '',
        lastTargetPdb: '',
        queryData: null,
        targetData: null,
    }),
    props: {
        alignment: { type: Object, required: true },
        lineLen: { type: Number, required: true },
        queryPdb: { type: String, required: true },
        targetPdb: { type: String, required: true },
        strRepr: { type: String, default: 'ribbon' },
        motifRepr: { type: String, default: 'licorice' },
        qmotifClr: { type: String, default: '#1E88E5' },
        tmotifClr: { type: String, default: '#FFC107' },
        qClr: { type: String, default: '#A5CFF5' },
        tClr: { type: String, default: '#FFE699' },
        qnullClr: { type: String, default: '#A5CFF5' },
        tnullClr: { type: String, default: '#FFE699' },
        bgColorLight: { type: String, default: 'white' },
        bgColorDark: { type: String, default: '#1E1E1E' },
        autoViewTime: { type: Number, default: 100 },
    },
    computed: {
        tmPanelBindings() {
            return this.isFullscreen ? { style: 'margin-top: 10px; font-size: 2em; line-height: 2em' } : {};
        },
    },
    methods: {
        handleToggleQuery() {
            if (!this.stage) return;
            if (__LOCAL__) {
                this.showQuery = this.showQuery === 0 ? 1 : 0;
            } else {
                this.showQuery = this.showQuery === 2 ? 0 : this.showQuery + 1;
            }
        },
        handleToggleTarget() {
            if (!this.stage) return;
            if (__LOCAL__) {
                this.showTarget = this.showTarget === 0 ? 1 : 0;
            } else {
                this.showTarget = this.showTarget === 2 ? 0 : this.showTarget + 1;
            }
        },
        handleResetView() {
            if (!this.stage) return;
            this.focusSelection();
        },
        async handleMakeImage() {
            if (!this.stage) return;
            const wasSpinning = this.isSpinning;
            this.isSpinning = false;
            const title = this.alignment?.target?.replace(/\.(pdb|cif)(\.gz)?$/i, '') || 'folddisco';
            const blob = await this.stage.makeImage();
            if (blob) {
                downloadBlob(blob, `folddisco-${title}.png`);
            }
            this.isSpinning = wasSpinning;
        },
        handleMakePDB() {
            const qPDB = extractAtomLines(this.lastQueryPdb);
            const tPDB = extractAtomLines(this.lastTargetPdb);
            if (qPDB.length === 0 && tPDB.length === 0) return;
            const title = `folddisco-${(this.alignment?.target || 'target').replace(/\.(pdb|cif)(\.gz)?$/i, '')}`;
            let result = null;
            if (qPDB.length > 0 && tPDB.length > 0) {
                result =
`TITLE     ${title}
REMARK     This file was generated by the Foldseek webserver:
REMARK       https://search.foldseek.com/folddisco
MODEL        1
${qPDB.join('\n')}
ENDMDL
MODEL        2
${tPDB.join('\n')}
ENDMDL
END
`;
            } else {
                result =
`TITLE     ${title}
REMARK     This file was generated by the Foldseek webserver:
REMARK       https://search.foldseek.com/folddisco
MODEL        1
${tPDB.join('\n')}
ENDMDL
END
`;
            }
            downloadBlob(new Blob([result], { type: 'text/plain' }), `${title}.pdb`);
        },
        focusSelection() {
            const target = this.queryStructureRef || this.targetStructureRef;
            if (!target) {
                this.stage.focusLoci(null, this.autoViewTime);
                return;
            }
            const expr = this.queryStructureRef
                ? (buildMotifExpression(this.alignment.queryresidues) || MS.struct.generator.all())
                : (buildMotifExpression(this.alignment.targetresidues) || MS.struct.generator.all());
            const loci = getSelectionLoci(expr, target);
            this.stage.focusLoci(loci, this.autoViewTime);
        },
        async loadStructureData() {
            const query = processPdb(this.queryPdb);
            const target = processPdb(this.targetPdb);
            if (!query || !target) {
                this.queryData = null;
                this.targetData = null;
                this.lastQueryPdb = '';
                this.lastTargetPdb = '';
                return;
            }

            this.queryData = { data: query.data, format: query.format, label: 'query' };
            this.targetData = { data: target.data, format: target.format, label: 'target' };
            this.lastQueryPdb = query.format === 'pdb' ? query.data : '';

            if (target.format === 'pdb' && this.alignment?.tmat && this.alignment?.umat) {
                const t = this.alignment.tmat.split(',').map(x => Number.parseFloat(x));
                const uFlat = this.alignment.umat.split(',').map(x => Number.parseFloat(x));
                if (t.length === 3 && uFlat.length === 9 && t.every(Number.isFinite) && uFlat.every(Number.isFinite)) {
                    const u = [
                        [uFlat[0], uFlat[1], uFlat[2]],
                        [uFlat[3], uFlat[4], uFlat[5]],
                        [uFlat[6], uFlat[7], uFlat[8]],
                    ];
                    const transformed = transformPdb(target.data, t, u);
                    this.targetData = { data: transformed, format: 'pdb', label: 'target' };
                    this.lastTargetPdb = transformed;
                    return;
                }
            }
            this.lastTargetPdb = target.format === 'pdb' ? target.data : '';
        },
        async renderStructures(focus = true) {
            if (!this.stage || !this.stageReady) {
                return;
            }
            const token = ++this.renderToken;
            await this.stageReady;
            if (token !== this.renderToken) {
                return;
            }

            await this.stage.clear();
            this.queryStructureRef = null;
            this.targetStructureRef = null;

            if (!this.queryData || !this.targetData) {
                return;
            }

            const backboneRepr = normalizeRepresentation(this.strRepr, 'cartoon');
            const motifRepr = normalizeRepresentation(this.motifRepr, 'ball-and-stick');

            const qMotif = buildMotifExpression(this.alignment.queryresidues);
            const tMotif = buildMotifExpression(this.alignment.targetresidues);
            const qChains = buildMotifChainExpression(this.alignment.queryresidues);
            const tChains = buildMotifChainExpression(this.alignment.targetresidues);

            const queryStructure = await this.stage.loadStructure(this.queryData);
            const targetStructure = await this.stage.loadStructure(this.targetData);
            this.queryStructureRef = queryStructure;
            this.targetStructureRef = targetStructure;

            const qMotifColor = toMolstarColor(this.qmotifClr, 0x1e88e5);
            const tMotifColor = toMolstarColor(this.tmotifClr, 0xffc107);
            const qChainColor = toMolstarColor(this.qClr, 0xa5cff5);
            const tChainColor = toMolstarColor(this.tClr, 0xffe699);
            const qOutsideColor = toMolstarColor(this.qnullClr, 0xa5cff5);
            const tOutsideColor = toMolstarColor(this.tnullClr, 0xffe699);

            if (this.showQuery >= 1 && qChains) {
                const queryChainComp = await this.stage.createComponentFromExpression(queryStructure, qChains, `query-chains-${token}`);
                if (queryChainComp) {
                    await this.stage.addRepresentation(queryChainComp, {
                        type: backboneRepr,
                        color: 'uniform',
                        colorParams: { value: qChainColor },
                        typeParams: { alpha: 0.8 },
                    });
                }
            }
            if (this.showQuery === 2 && qChains) {
                const outsideExpr = MS.struct.modifier.exceptBy({ 0: MS.struct.generator.all(), by: qChains });
                const queryOutsideComp = await this.stage.createComponentFromExpression(queryStructure, outsideExpr, `query-outside-${token}`);
                if (queryOutsideComp) {
                    await this.stage.addRepresentation(queryOutsideComp, {
                        type: backboneRepr,
                        color: 'uniform',
                        colorParams: { value: qOutsideColor },
                        typeParams: { alpha: 0.35 },
                    });
                }
            }
            if (qMotif) {
                const queryMotifComp = await this.stage.createComponentFromExpression(queryStructure, qMotif, `query-motif-${token}`);
                if (queryMotifComp) {
                    await this.stage.addRepresentation(queryMotifComp, {
                        type: motifRepr,
                        color: 'uniform',
                        colorParams: { value: qMotifColor },
                    });
                }
            }

            if (this.showTarget >= 1 && tChains) {
                const targetChainComp = await this.stage.createComponentFromExpression(targetStructure, tChains, `target-chains-${token}`);
                if (targetChainComp) {
                    await this.stage.addRepresentation(targetChainComp, {
                        type: backboneRepr,
                        color: 'uniform',
                        colorParams: { value: tChainColor },
                        typeParams: { alpha: 0.8 },
                    });
                }
            }
            if (this.showTarget === 2 && tChains) {
                const outsideExpr = MS.struct.modifier.exceptBy({ 0: MS.struct.generator.all(), by: tChains });
                const targetOutsideComp = await this.stage.createComponentFromExpression(targetStructure, outsideExpr, `target-outside-${token}`);
                if (targetOutsideComp) {
                    await this.stage.addRepresentation(targetOutsideComp, {
                        type: backboneRepr,
                        color: 'uniform',
                        colorParams: { value: tOutsideColor },
                        typeParams: { alpha: 0.35 },
                    });
                }
            }
            if (tMotif) {
                const targetMotifComp = await this.stage.createComponentFromExpression(targetStructure, tMotif, `target-motif-${token}`);
                if (targetMotifComp) {
                    await this.stage.addRepresentation(targetMotifComp, {
                        type: motifRepr,
                        color: 'uniform',
                        colorParams: { value: tMotifColor },
                    });
                }
            }

            if (focus) {
                this.focusSelection();
            }
        },
    },
    watch: {
        showQuery() {
            this.renderStructures(false);
        },
        showTarget() {
            this.renderStructures(false);
        },
        alignment: {
            deep: true,
            handler() {
                this.loadStructureData().then(() => this.renderStructures(true));
            },
        },
        queryPdb() {
            this.loadStructureData().then(() => this.renderStructures(true));
        },
        targetPdb() {
            this.loadStructureData().then(() => this.renderStructures(true));
        },
    },
    async mounted() {
        if (!this.alignment || !this.queryPdb || !this.targetPdb) {
            return;
        }
        await this.stageReady;
        await this.loadStructureData();
        await this.renderStructures(true);
    },
};
</script>

<style scoped>
.structure-panel {
    width: 100%;
    height: 100%;
    position: relative;
}
.structure-wrapper {
    width: 500px;
    height: 400px;
    margin: 0 auto;
}

.structure-viewer {
    width: 100%;
    height: 100%;
    overflow: hidden;
    padding: 0;
}

.structure-viewer canvas {
    border-radius: 2px;
}

@media screen and (max-width: 960px) {
    .structure-panel {
        display: flex;
        flex-direction: column-reverse;
    }
    .structure-wrapper {
        align-self: center;
        padding-bottom: 1em;
    }
}

@media screen and (min-width: 961px) {
    .structure-wrapper {
        padding-left: 2em;
    }
}
</style>
