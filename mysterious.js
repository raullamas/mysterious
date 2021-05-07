var settings = getSettings(undefined, 15);
var errors = getErrors();
var protoSpecim = getProtoSpecim();

console.log(settings);

//

function findMostRelated (specims) {
    if (!Array.isArray(specims) || specims.length < 2){
        return 'Error: ' + errors.isNotArrayTwoEl;
    } else if (specims.some(isInvalidSpecim)) {
        return 'Error: ' + errors.invalidSpecims;
    }
    return specims.reduce(compareToAllOther, {'top relatedness': 0, IDs: [], strands: []});

    function compareToAllOther (mostRelated, specim) {
        let toCompareTo = specims.slice(specims.indexOf(specim) + 1);
        for (let i = 0; i < toCompareTo.length; i++) {
            let compared = specim.compare(toCompareTo[i])
            if (compared.relatedness > mostRelated['top relatedness']) {
                mostRelated.IDs = [compared.IDs];
                mostRelated.strands = [compared.strands];
                mostRelated['top relatedness'] = compared.relatedness;
            } else if (compared.relatedness === mostRelated['top relatedness']) {
                mostRelated.IDs.push(compared.IDs);
                mostRelated.strands.push(compared.strands);
            }
        }
        return mostRelated;
    };
}

function getErrors () {
    return {
        invalidID: `ID must be a numeric value between 0 and ${settings.poolSize - 1} (inclusive).`,
        invalidStrand: `${settings.nucAcid === 'DNA' ? 'DNA' : 'RNA'} strand must be composed exclusively of bases ${settings.nucAcid === 'DNA' ? ['A', 'T', 'C', 'G'] : ['A', 'U', 'C', 'G']} and must be ${settings.strandLength} bases long.`,
        poolIsFull: `There are already ${settings.poolSize} unique specimens on record. Update pool size to create new specimens.\n\n(Maximum number specimens is determined by raising number of nucleobases (4) to the power of strand length (${settings.strandLength}).)`,
        duplicateID (id) {
            return `ID ${id} is already on record. Choose a different ID.`;
        },
        duplicateStrand (strand) {
            return `${settings.nucAcid === 'DNA' ? 'DNA' : 'RNA'} strand ${strand} is already on record. Choose a different strand.`;
        },
        invalidSpecim: 'Specimen to compare to must be created through command "newSpecim".',
        invalidNumMutations: `Number of mutations must be a numeric value between 1 and ${settings.strandLength} (inclusive).`,
        invalidSturdinessPcent: 'Sturdiness percentage must be a numeric value between 0 and 100 (inclusive).',
        invalidNumSturdy () {
            return `Number of specimens must be a numeric value between 1 and ${settings.poolSize - settings.IDs.size} (number of available specimens as per pool size), inclusive.`
        },
        isNotArrayTwoEl: 'Specimens to evaluate must be passed inside an array of at least 2 elements.',
        invalidSpecims: 'All specimens to evaluate must be created through command "newSpecim".'
    }
}

function getNRandSturdy(num = 30, pcent = 60) {
    if (typeof num !== 'number' || num < 1 || num > settings.poolSize - settings.IDs.size) {
        return 'Error: ' + errors.invalidNumSturdy();
    } else if (typeof pcent !== 'number' || pcent < 0 || pcent > 100) {
        return 'Error: ' + errors.invalidSturdinessPcent;
    }
    var sturdy = [];
    while (sturdy.length < num) {
        let randSpecim = newSpecim();
        if (randSpecim.isSturdy(pcent)) {
            sturdy.push(randSpecim);
        } else {
            settings.IDs.delete(randSpecim.ID);
            settings.strands.delete(randSpecim.strand);
        }
    }
    return sturdy;
}

function getProtoSpecim() {
    return {
        compare(specim) {
            if (isInvalidSpecim(specim)) {
                return 'Error: ' + errors.invalidSpecim;
            }
            let simil = 0;
            for (let i = 0, j = 0; i < this.strand.length; i++, j++) {
                if (this.strand[i] === specim.strand[j]) {
                    simil++;
                }
            };
            return { 
                IDs: [this.ID, specim.ID],
                strands: [this.strand, specim.strand],
                relatedness: +(simil / this.strand.length * 100).toFixed(2)
            }
        },
        complement() {
            return Array.from(this.strand).map(complemBase);

            function complemBase(base) {
                switch (base) {
                    case 'A':
                        return settings.nucAcid === 'DNA' ? 'T' : 'U';
                    case 'T':
                    case 'U':
                        return 'A';
                    case 'C':
                        return 'G';
                    case 'G':
                        return 'C';
                    default:
                        return 'x';
                }
            }
        },
        mutate(num = 5) {
            if (typeof num !== 'number' || num > settings.strandLength || num < 1) {
                return 'Error: ' + errors.invalidNumMutations;
            }
            let [newStrand, numSwaps, idxsOfSwaps] = [this.strand, 0, new Set()];
            while (numSwaps < num) {
                let [swapIdx, swappedBase] = [Math.floor(Math.random() * newStrand.length), randBase()];
                if (newStrand[swapIdx] === swappedBase || idxsOfSwaps.has(swapIdx)) {
                    continue;
                }
                newStrand = swapOneBase(newStrand, swapIdx, swappedBase);
                numSwaps++;
                idxsOfSwaps.add(swapIdx);
            }
            return newStrand;

            function swapOneBase(strand, idx, replacement) {
                return strand.substring(0, idx) + replacement + strand.substring(idx + 1);
            }
        }, 
        isSturdy(pcent = 60) {
            if (typeof pcent !== 'number' || pcent > 100 || pcent < 0) {
                return 'Error: ' + errors.invalidSturdinessPcent;
            }
            return Array.from(this.strand).reduce(incrementIf, 0) >= pcent/100 * settings.strandLength;

            function incrementIf(acc, base) {
                if (base === 'C' || base === 'G') {
                    acc++;
                }
                return acc;
            }
        }
    }
}

function getSettings (nucAcid = 'DNA', strandLength = 15) {
    if (nucAcid !== 'DNA' && nucAcid !== 'RNA') {
        return 'Error: ' + 'Nucleic acid must be either \'DNA\' or \'RNA\'.';
    } else if (typeof strandLength !== 'number' || strandLength > 15 || strandLength < 2) {
        return 'Error: ' + 'Strand length must be a numeric value between 2 and 15 (inclusive)';
    }
    return {
        nucAcid,
        get bases () {
            return this.nucAcid === 'DNA' ? ['A', 'T', 'C', 'G'] : ['A', 'U', 'C', 'G']
        },
        strandLength,
        get poolSize () {
            return Math.pow(4, this.strandLength);
        },
        IDs: new Set(),
        strands: new Set(),
    }
}

function isInvalidSpecim (specim) {
    return Object.getPrototypeOf(Object(specim)) !== protoSpecim;
}

function newSpecim(ID = randID(), strand = randStrand()) {
    if (typeof ID !== 'number' || ID >= settings.poolSize || ID < 0) {
        return 'Error: ' + errors.invalidID;
    } else if (typeof strand !== 'string' || strand.length !== settings.strandLength || isInvalidStrand(strand)) {
        return 'Error: ' + errors.invalidStrand;
    };
    if (settings.IDs.size >= settings.poolSize || settings.strands.size >= settings.poolSize) {
        return 'Error: ' + errors.poolIsFull;
    };
    if (settings.IDs.has(ID)) {
        return 'Error: ' + errors.duplicateID(ID);
    } else if (settings.strands.has(strand)) {
        return 'Error: ' + errors.duplicateStrand(strand);
    };
    settings.IDs.add(ID), settings.strands.add(strand);
    let specim = Object.create(protoSpecim);
    specim.ID = ID, specim.strand = strand;
    return specim;

    function isInvalidStrand(strand) {
        let strandArr = Array.from(strand);
        return !strandArr.every(base => settings.bases.includes(base));
    }
}

function randBase() {
    const bases = settings.bases;
    return bases[Math.floor(Math.random() * 4)];
}

function randID() {
    if (settings.IDs.size < settings.poolSize) {
        let ID = Math.floor(Math.random() * settings.poolSize);
        if (settings.IDs.has(ID)) {
            return randID();
        };
        return ID;
    }
}

function randStrand() {
    if (settings.strands.size < settings.poolSize) {
        let newStrand = '';
        for (let i = 0; i < settings.strandLength; i++) {
            newStrand += randBase();
        }
        if (settings.strands.has(newStrand)) {
            return randStrand();
        };
        return newStrand;
    }
}