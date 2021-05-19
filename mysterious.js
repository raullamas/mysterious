var settings = getSettings(' rnA ');
var errors = getErrors();
var protoSpecim = getProtoSpecim();

console.log(settings);
console.log(findMostRelated([newSpecim(), {}]));
console.log(settings);

//

function findMostRelated (specims) {
    if (!Array.isArray(specims)) {
        return 'Error: ' + errors.notArray('Specimens to evaluate');
    }
    specims = specims.flat();
    if (specims.length < 2) {
        return 'Error: ' + errors.needMoreEl('Specimens to evaluate', 2);
    } else if (specims.some(isInvalidSpecim)) {
        return 'Error: ' + errors.invalidSpecims('Specimens to evaluate');
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
        notString (arg) {
            return `${arg} must be passed as a string: Use quotation marks (\'\').`;
        },
        notArray (arg) {
            return `${arg} must be passed as an array: Use brackets (\'[ ]\') and commas (\',\').`;
        },
        notNumber (arg) {
            return `${arg} must be a valid integer: Use only digits (0 - 9) and pass at least one digit.`;
        },
        offRange (arg, lowLim, upLim) {
            return `${arg} must be greater than ${lowLim - 1} and less than ${upLim + 1}: Update value accordingly.`;
        },
        notStrand (arg, bases) {
            return `${arg} must consist only of bases ${bases.join(', ')}: Update strand accordingly.`;
        },
        invalidLength (arg, length) {
            return `${arg} must be ${length} characters long: Update strand accordingly.`;
        },
        needMoreEl (arg, numEl) {
            return `${arg} must be an array of at least ${numEl} elements: Update array accordingly.`;
        },
        poolFull (arg, upLim) {
            return `${arg} is already maxed out at ${upLim}: Increase strand length to create new specimens.`;
        },
        duplicated (arg, datum) {
            return `${arg} ${datum} is already on record: Choose a different ${arg.split(' ')[1]}.`;
        },
        invalidSpecim (arg) {
            return `${arg} must be created through command 'newSpecim': Update specimen accordingly.`;
        },
        invalidSpecims (arg) {
            return `${arg} must all be created through command 'newSpecim': Update specimens accordingly.`;
        }
    }
}

function getNRandSturdy(num = '30', pcent = '60') {
    if (typeof num !== 'string') {
        return 'Error: ' + errors.notString('Number of specimens');
    } else if (isInvalidNum(num.trim())) {
        return 'Error: ' + errors.notNumber('Number of specimens');
    };
    num = parseInt(num.trim());
    if (num > settings.poolSize - settings.IDs.size || num < 1) {
        return 'Error: ' + errors.offRange('Number of specimens', 1, settings.poolSize - settings.IDs.size);
    };
    
    if (typeof pcent !== 'string') {
        return 'Error: ' + errors.notString('Sturdiness percentage');
    } else if (isInvalidNum(pcent.trim())) {
        return 'Error: ' + errors.notNumber('Sturdiness percentage');
    };
    pcent = parseInt(pcent.trim());
    if (pcent > 100 || pcent < 0) {
        return 'Error: ' + errors.offRange('Sturdiness percentage', 0, 100);
    };

    var sturdy = [];
    while (sturdy.length < num) {
        let randSpecim = newSpecim();
        if (randSpecim.isSturdy('' + pcent)) {
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
                return 'Error: ' + errors.invalidSpecim('Specimen to compare to');
            };

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
        mutate(num = '5') {
            if (typeof num !== 'string') {
                return 'Error: ' + errors.notString('Number of mutations');
            } else if (isInvalidNum(num.trim())) {
                return 'Error: ' + errors.notNumber('Number of mutations');
            };
            num = parseInt(num.trim());
            if (num > settings.strandLength || num < 1) {
                return 'Error: ' + errors.offRange('Number of mutations', 1, settings.strandLength);
            };

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
        isSturdy(pcent = '60') {
            if (typeof pcent !== 'string') {
                return 'Error: ' + errors.notString('Sturdiness percentage');
            } else if (isInvalidNum(pcent.trim())) {
                return 'Error: ' + errors.notNumber('Sturdiness percentage');
            };
            pcent = parseInt(pcent.trim());
            if (pcent > 100 || pcent < 0) {
                return 'Error: ' + errors.offRange('Sturdiness percentage', 0, 100);
            }

            return Array.from(this.strand).reduce(incrementIf, 0) >= parseInt(pcent)/100 * settings.strandLength;

            function incrementIf(acc, base) {
                if (base === 'C' || base === 'G') {
                    acc++;
                }
                return acc;
            }
        }
    }
}

function getSettings (nucAcid = 'DNA', strandLength = '15') {
    if (typeof nucAcid !== 'string') {
        return 'Error: ' + `${'Nucleic acid'} must be passed as a string: Please surround with quotation marks (\'\').`;
    };
    nucAcid = nucAcid.trim().toUpperCase();
    if (nucAcid !== 'DNA' && nucAcid !== 'RNA') {
        return 'Error: ' + `${'Nucleic acid'} must be either ${'DNA'} or ${'RNA'}.`;
    }
    if (typeof strandLength !== 'string') {
        return 'Error: ' + `${'Strand length'} must be passed as a string: Please surround with quotation marks (\'\').`;
    } else if (isInvalidNum(strandLength.trim())) {
        return 'Error: ' + `${'Strand length'} must be a valid integer value: Remove non-numeric characters and/or fractional components from input, and make sure it has at least 1 digit.`;
    };
    strandLength = parseInt(strandLength.trim());
    if (strandLength > 15 || strandLength < 2) {
        return 'Error: ' + `${'Strand length'} must be greater than than ${2 - 1} and less than ${15 + 1}.`;
    };
    
    return {
        nucAcid,
        get bases () {
            return this.nucAcid === 'DNA' ? ['A', 'T', 'C', 'G'] : ['A', 'U', 'C', 'G'];
        },
        strandLength,
        get poolSize () {
            return Math.pow(4, this.strandLength);
        },
        IDs: new Set(),
        strands: new Set(),
    }
}

function isInvalidNum (num) {
    return Array.from(num).some(char => Object.is(parseInt(char), NaN)) || num === '';
}

function isInvalidSpecim (specim) {
    return Object.getPrototypeOf(Object(specim)) !== protoSpecim;
}

function newSpecim(ID = randID(), strand = randStrand()) {
    if (settings.IDs.size >= settings.poolSize) {
        return 'Error: ' + errors.poolFull('Number of specimens on record', settings.poolSize);
    } else if (typeof ID !== 'string') {
        return 'Error: ' + errors.notString('Specimen ID');
    } else if (isInvalidNum(ID.trim())) {
        return 'Error: ' + errors.notNumber('Specimen ID');
    };
    ID = pad(ID.trim());
    if (parseInt(ID) > settings.poolSize || parseInt(ID) < 1) {
        return 'Error: ' + errors.offRange('Specimen ID', 1, settings.poolSize);
    } else if (settings.IDs.has(ID)) {
        return 'Error: ' + errors.duplicated('Specimen ID', ID);
    };
    if (typeof strand !== 'string') {
        return 'Error: ' + errors.notString(`${settings.nucAcid} strand`);
    };
    strand = strand.trim().toUpperCase();
    if (isInvalidStrand(strand)) {
        return 'Error: ' + errors.notStrand(`${settings.nucAcid} strand`, settings.bases);
    } else if (strand.length !== settings.strandLength) {
        return 'Error: ' + errors.invalidLength(`${settings.nucAcid} strand`, settings.strandLength);
    } else if (settings.strands.has(strand)) {
        return 'Error: ' + errors.duplicated(`${settings.nucAcid} strand`, strand);
    };

    settings.IDs.add(ID), settings.strands.add(strand);
    let specim = Object.create(protoSpecim);
    specim.ID = ID, specim.strand = strand;
    return specim;

    function isInvalidStrand(strand) {
        return !Array.from(strand).every(base => settings.bases.includes(base));
    }
}

function pad (ID) {
    return ID.padStart(settings.poolSize.toString().length, '0');
}

function randBase() {
    return settings.bases[Math.floor(Math.random() * 4)];
}

function randID() {
    if (settings.IDs.size < settings.poolSize) {
        let ID = pad('' + (Math.floor(Math.random() * settings.poolSize) + 1));
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