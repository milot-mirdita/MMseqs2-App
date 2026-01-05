import { PluginContext } from 'molstar/lib/mol-plugin/context';
import { Color } from 'molstar/lib/mol-util/color';

const DEFAULT_SPIN_SPEED = 0.05;

function parseHexColor(value, fallback = 0xffffff) {
    if (typeof value === 'number' && Number.isFinite(value)) {
        return Color(value);
    }
    if (typeof value === 'string') {
        const trimmed = value.trim();
        if (trimmed.length === 0) {
            return Color(fallback);
        }
        const hex = trimmed.startsWith('#') ? trimmed.slice(1) : trimmed;
        const parsed = Number.parseInt(hex, 16);
        if (Number.isFinite(parsed)) {
            return Color(parsed);
        }
    }
    return Color(fallback);
}

export class MolstarStage {
    constructor(container, params = {}) {
        this.container = container;
        this.params = params;
        this.plugin = null;
        this.canvas = null;
        this.spinSpeed = DEFAULT_SPIN_SPEED;
    }

    async init() {
        if (!this.container) return;
        this.plugin = new PluginContext({ actions: [], behaviors: [], animations: [], config: [] });
        await this.plugin.init();

        if (this.canvas && this.canvas.parentElement) {
            this.canvas.parentElement.removeChild(this.canvas);
        }
        this.canvas = document.createElement('canvas');
        this.canvas.setAttribute('aria-hidden', 'true');
        this.container.appendChild(this.canvas);

        const ok = await this.plugin.initViewerAsync(this.canvas, this.container);
        if (ok === false) {
            throw new Error('Mol* viewer initialization failed');
        }

        this.plugin.canvas3d?.setProps({
            transparentBackground: true,
            sceneRadiusFactor: 3.0,
            camera: {
                helper: { axes: { name: 'off', params: {} } },
                fov: 60,
            },
            cameraClipping: {
                radius: 0,
                far: false,
                minNear: -1000,
            },
        });

        if (this.params.backgroundColor) {
            this.setBackground(this.params.backgroundColor);
        }

        this.plugin.canvas3d?.handleResize();
    }

    async clear() {
        if (!this.plugin) return;
        await this.plugin.clear();
    }

    async loadStructure({ data, format, label }) {
        if (!this.plugin) return null;
        const raw = await this.plugin.builders.data.rawData({ data, label }, { state: { isGhost: true } });
        // console.log(raw)
        const trajectory = await this.plugin.builders.structure.parseTrajectory(raw, format);
        const model = await this.plugin.builders.structure.createModel(trajectory);
        return this.plugin.builders.structure.createStructure(model);
    }

    async createComponentFromExpression(structure, expression, key) {
        if (!this.plugin || !structure || !expression) return null;
        return this.plugin.builders.structure.tryCreateComponentFromExpression(structure, expression, key);
    }

    async createComponentStatic(structure, type = 'all') {
        if (!this.plugin || !structure) return null;
        return this.plugin.builders.structure.tryCreateComponentStatic(structure, type);
    }

    async addRepresentation(component, props, options) {
        if (!this.plugin || !component) return null;
        return this.plugin.builders.structure.representation.addRepresentation(component, props, options);
    }

    async addDistance(lociA, lociB, options = {}) {
        if (!this.plugin || !lociA || !lociB) return null;
        return this.plugin.managers.structure.measurement.addDistance(lociA, lociB, options);
    }

    async remove(ref) {
        if (!this.plugin || !ref) return;
        const targetRef = typeof ref === 'string' ? ref : ref.ref;
        if (!targetRef) return;
        const update = this.plugin.state.data.build().delete(targetRef);
        if (update.editInfo.count === 0) return;
        await update.commit();
    }

    onClick(handler) {
        if (!this.plugin?.behaviors?.interaction?.click) return () => {};
        const sub = this.plugin.behaviors.interaction.click.subscribe(handler);
        return () => sub.unsubscribe();
    }

    setBackground(color) {
        if (!this.plugin?.canvas3d) return;
        const renderer = this.plugin.canvas3d.props.renderer;
        this.plugin.canvas3d.setProps({
            renderer: {
                ...renderer,
                backgroundColor: parseHexColor(color),
            },
        });
    }

    setSpin(enabled) {
        if (!this.plugin?.canvas3d) return;
        const trackball = this.plugin.canvas3d.props.trackball;
        this.plugin.canvas3d.setProps({
            trackball: {
                ...trackball,
                animate: enabled
                    ? { name: 'spin', params: { speed: this.spinSpeed } }
                    : { name: 'off', params: {} },
            },
        });
    }

    focusLoci(loci, durationMs = 0) {
        if (!this.plugin) return;
        if (loci) {
            this.plugin.managers.camera.focusLoci(loci, { durationMs });
        } else {
            this.plugin.managers.camera.reset();
        }
    }

    handleResize() {
        this.plugin?.canvas3d?.handleResize();
    }

    toggleFullscreen(element) {
        const target = element || this.container;
        if (!target) return;
        if (!document.fullscreenElement && target.requestFullscreen) {
            target.requestFullscreen();
        } else if (document.fullscreenElement && document.exitFullscreen) {
            document.exitFullscreen();
        }
    }

    async makeImage() {
        if (!this.canvas) return null;
        return new Promise(resolve => {
            this.canvas.toBlob(blob => resolve(blob), 'image/png');
        });
    }

    dispose() {
        if (this.plugin) {
            this.plugin.dispose();
            this.plugin = null;
        }
        if (this.canvas && this.canvas.parentElement) {
            this.canvas.parentElement.removeChild(this.canvas);
        }
        this.canvas = null;
    }
}
