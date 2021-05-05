var onRecord = getOnRecord();
var errors = getErrors();
var projectMethods = coreMethods();

let manySturdy = getNRandSturdy(5, 100);
console.log(manySturdy);
console.log(onRecord);

//
function coreMethods() {
    return {
        compareDna(specim) {
            if (isInvalidSpecim(specim)) {
                return 'Error: ' + errors.invalidSpecim;
            }
            let simil = 0;
            for (let i = 0, j = 0; i < this.dnaSeq.length; i++, j++) {
                if (this.dnaSeq[i] === specim.dnaSeq[j]) {
                    simil++;
                }
            };
            return { 
                iDs: [this.specimId, specim.specimId],
                dnaSeqs: [this.dnaSeq, specim.dnaSeq],
                compatibility: +(simil / this.dnaSeq.length * 100).toFixed(2)
            }
        },
        complemDnaSeq() {
            return Array.from(this.dnaSeq).map(complemDnaBase);

            function complemDnaBase(base) {
                switch (base) {
                    case 'A':
                        return onRecord.nucAcid === 'DNA' ? 'T' : 'U';
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
        mutate(num = 1) {
            if (typeof num !== 'number' || num > this.dnaSeq.length || num < 1) {
                return 'Error: ' + errors.invalidNumMutations;
            }
            let [mutatedSeq, numSwaps, idxsOfSwaps] = [this.dnaSeq, 0, new Set()];
            while (numSwaps < num) {
                let [swapIdx, swappedBase] = [Math.floor(Math.random() * mutatedSeq.length), randDnaBase()];
                if (mutatedSeq[swapIdx] === swappedBase || idxsOfSwaps.has(swapIdx)) {
                    continue;
                }
                mutatedSeq = swapOneChar(mutatedSeq, swapIdx, swappedBase);
                numSwaps++;
                idxsOfSwaps.add(swapIdx);
            }
            return mutatedSeq;

            function swapOneChar(str, idx, replacement) {
                return str.substring(0, idx) + replacement + str.substring(idx + 1);
            }
        }, 
        isSturdy(pcent = 60) {
            if (typeof pcent !== 'number' || pcent > 100 || pcent < 0) {
                return 'Error: ' + errors.invalidSturdinessPcent;
            }
            return Array.from(this.dnaSeq).reduce(reducer, 0) >= pcent/100 * onRecord.strandLength;

            function reducer(acc, curVal) {
                if (curVal === 'C' || curVal === 'G') {
                    acc++;
                }
                return acc;
            }
        }
    }
};

function findMostRelated (specimens) {
    if (!Array.isArray(specimens) || specimens.length < 2){
        return 'Error: ' + errors.isNotArray;
    } else if (specimens.some(isInvalidSpecim)) {
        return 'Error: ' + errors.invalidSpecims;
    }
    return specimens.reduce(reducer1, {'top compatibility': 0, iDs: [], dnaSeqs: []});

    function reducer1 (acc, curVal) {
        let notYetCompared = specimens.slice(specimens.indexOf(curVal));
        for (let i = 1; i < notYetCompared.length; i++) {
            let compared = curVal.compareDna(notYetCompared[i])
            if (compared.compatibility > acc['top compatibility']) {
                acc.iDs = [compared.iDs];
                acc.dnaSeqs = [compared.dnaSeqs];
                acc['top compatibility'] = compared.compatibility;
            } else if (compared.compatibility === acc['top compatibility']) {
                acc.iDs.push(compared.iDs);
                acc.dnaSeqs.push(compared.dnaSeqs);
            }
        }
        return acc;
    };
}

function getErrors () {
    return {
        invalidID: `Specimen ID must be a numeric value between 0 and ${onRecord.poolSize - 1} (inclusive).`,
        invalidDnaSeq: `Strand must be composed exclusively of bases ${onRecord.nucAcid === 'DNA' ? ['A', 'T', 'C', 'G'] : ['A', 'U', 'C', 'G']} and must be ${onRecord.strandLength} bases long.`,
        poolIsFull: `There are already ${onRecord.poolSize} unique specimens on record. Update pool size to create new specimens.\n\n(Maximum number specimens is determined by raising number of nucleobases (4) to the power of user-defined DNA sequence length (${onRecord.strandLength}).)`,
        duplicateID (id) {
            return `Specim ID ${id} is already on record. Choose a different ID.`;
        },
        duplicateDNASeq (dnaSeq) {
            return `DNA sequence ${dnaSeq} is already on record. Choose a different DNA sequence.`;
        },
        invalidSpecim: 'Specimen to compare to must be created through command "newSpecim".',
        invalidNumMutations: `Number of mutations must be a numeric value between 1 and ${onRecord.strandLength} (inclusive).`,
        invalidSturdinessPcent: 'Sturdiness percentage must be a numeric value between 0 and 100 (inclusive).',
        invalidNumRandSturdy () {
            return `Number of specimens must be a numeric value between 1 and ${onRecord.poolSize - onRecord.iDs.size} (number of available specimens as per pool size), inclusive.`
        },
        isNotArray: 'Specimens to evaluate must be passed as an array of at least 2 elements.',
        invalidSpecims: 'All specimens to evaluate must be created through command "newSpecim".'
    }
}

function getNRandSturdy(num, pcent = 60) {
    if (typeof num !== 'number' || num < 1 || num > onRecord.poolSize - onRecord.iDs.size) {
        return 'Error: ' + errors.invalidNumRandSturdy();
    } else if (typeof pcent !== 'number' || pcent < 0 || pcent > 100) {
        return 'Error: ' + errors.invalidSturdinessPcent;
    }
    var sturdy = [];
    while (sturdy.length < num) {
        let randSpecim = newSpecim();
        if (randSpecim.isSturdy(pcent)) {
            sturdy.push(randSpecim);
        } else {
            onRecord.iDs.delete(randSpecim.specimId);
            onRecord.dnaSeqs.delete(randSpecim.dnaSeq);
        }
    }
    return sturdy;
}

function getOnRecord (nucAcid = 'DNA', strandLength = 15) {
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
        iDs: new Set(),
        dnaSeqs: new Set(),
    }
}

function isInvalidSpecim (specim) {
    return Object.getPrototypeOf(Object(specim)) !== projectMethods;
}

function newSpecim(specimId = randSpecimId(), dnaSeq = randDnaSeq()) {
    if (typeof specimId !== 'number' || specimId >= onRecord.poolSize || specimId < 0) {
        return 'Error: ' + errors.invalidID;
    } else if (typeof dnaSeq !== 'string' || dnaSeq.length !== onRecord.strandLength || isInvalidSeq(dnaSeq)) {
        return 'Error: ' + errors.invalidDnaSeq;
    };
    if (onRecord.iDs.size >= onRecord.poolSize || onRecord.dnaSeqs.size >= onRecord.poolSize) {
        return 'Error: ' + errors.poolIsFull;
    };
    if (onRecord.iDs.has(specimId)) {
        return 'Error: ' + errors.duplicateID(specimId);
    } else if (onRecord.dnaSeqs.has(dnaSeq)) {
        return 'Error: ' + errors.duplicateDNASeq(dnaSeq);
    };
    onRecord.iDs.add(specimId), onRecord.dnaSeqs.add(dnaSeq);
    let specimen = Object.create(projectMethods);
    specimen.specimId = specimId, specimen.dnaSeq = dnaSeq;
    return specimen;

    function isInvalidSeq(str) {
        let arr = Array.from(str);
        return !arr.every(char => onRecord.bases.includes(char));
    }
}

function randDnaSeq() {
    if (onRecord.dnaSeqs.size < onRecord.poolSize) {
        let newStrand = '';
        for (let i = 0; i < onRecord.strandLength; i++) {
            newStrand += randDnaBase();
        }
        if (onRecord.dnaSeqs.has(newStrand)) {
            return randDnaSeq();
        };
        return newStrand;
    }
}

function randSpecimId() {
    if (onRecord.iDs.size < onRecord.poolSize) {
        let specimId = Math.floor(Math.random() * onRecord.poolSize);
        if (onRecord.iDs.has(specimId)) {
            return randSpecimId();
        };
        return specimId;
    }
}

function randDnaBase() {
    const bases = onRecord.bases;
    return bases[Math.floor(Math.random() * bases.length)];
}