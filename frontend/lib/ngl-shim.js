const unavailable = (name) => () => {
    throw new Error(`NGL stub: ${name} is not available`);
};

export class Stage {
    constructor() {
        unavailable('Stage')();
    }
}

export class Shape {
    constructor() {
        unavailable('Shape')();
    }
}

export class Selection {
    constructor() {
        unavailable('Selection')();
    }
}

export class Matrix4 {
    constructor() {
        unavailable('Matrix4')();
    }
}

export class Quaternion {
    constructor() {
        unavailable('Quaternion')();
    }
}

export class Vector3 {
    constructor() {
        unavailable('Vector3')();
    }
}

export class PdbWriter {
    constructor() {
        unavailable('PdbWriter')();
    }
    getData() {
        unavailable('PdbWriter.getData')();
    }
}

export const ColormakerRegistry = {
    addScheme: unavailable('ColormakerRegistry.addScheme'),
    addSelectionScheme: unavailable('ColormakerRegistry.addSelectionScheme'),
    hasScheme: unavailable('ColormakerRegistry.hasScheme'),
    removeScheme: unavailable('ColormakerRegistry.removeScheme'),
};

export const Color = (value) => value;

export const concatStructures = unavailable('concatStructures');
export const download = unavailable('download');
export const autoLoad = unavailable('autoLoad');
export const StructureComponent = unavailable('StructureComponent');
