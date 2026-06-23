// ── Init ──
document.addEventListener('DOMContentLoaded', () => {
    const kf = document.getElementById('ketcher-frame');
    if (kf) {
        kf.addEventListener('load', () => {
            setTimeout(() => window.scrollTo(0, 0), 50);
        });
    }
});

// ── Toast Notifications ──
function showToast(message, type, duration) {
    type = type || 'error';
    duration = duration || 4000;
    var container = document.getElementById('toast-container');
    var toast = document.createElement('div');
    toast.className = 'toast ' + type;
    toast.textContent = message;
    container.appendChild(toast);
    setTimeout(function () {
        toast.style.animation = 'toastOut 0.25s ease forwards';
        setTimeout(function () { toast.remove(); }, 250);
    }, duration);
}

// ── Tab Switching ──
function switchTab(tabId) {
    document.querySelectorAll('.tab-btn').forEach(function (btn) { btn.classList.remove('active'); });
    document.querySelectorAll('.tab-content').forEach(function (c) { c.classList.remove('active'); });
    document.querySelector('button[onclick="switchTab(\'' + tabId + '\')"]').classList.add('active');
    document.getElementById('tab-' + tabId).classList.add('active');
}

// ── Ketcher Transfer ──
async function getSmiles() {
    var kf = document.getElementById('ketcher-frame');
    if (kf && kf.contentWindow && kf.contentWindow.ketcher) {
        try {
            var smiles = await kf.contentWindow.ketcher.getSmiles();
            if (smiles) {
                document.getElementById('smiles-input').value = smiles;
                showToast('SMILES transferred', 'success', 2000);
            } else {
                showToast('Draw a molecule first.', 'info');
            }
        } catch (e) {
            showToast('Could not read from editor.');
        }
    } else {
        showToast('Editor is still loading.', 'info');
    }
}

// ── PubChem Lookup + Auto Predict ──
async function lookupAndPredict() {
    var cidInput = document.getElementById('pubchem-cid');
    var cid = cidInput.value.trim();
    if (!cid) {
        showToast('Enter a PubChem CID.', 'info');
        return;
    }

    var resultDiv = document.getElementById('pubchem-result');
    resultDiv.style.display = 'none';

    try {
        var response = await fetch('/herb-pred/pubchem/' + cid);
        var data = await response.json();

        if (!response.ok) {
            throw new Error(data.error || 'Lookup failed');
        }

        // Show compound info
        document.getElementById('pubchem-name').textContent = data.name || 'Unknown compound';
        document.getElementById('pubchem-formula').textContent = data.formula ? 'Formula: ' + data.formula : '';
        document.getElementById('pubchem-smiles').textContent = data.smiles;
        resultDiv.style.display = 'block';

        // Fill the SMILES box
        document.getElementById('smiles-input').value = data.smiles;

        // Auto-run prediction
        showToast('Fetched CID ' + cid + ' — running prediction…', 'success', 2000);
        await predictSingle();

    } catch (error) {
        showToast('PubChem error: ' + error.message);
    }
}

// ── Chart Instances ──
var radarChartInstance = null;
var fingerprintDonutChart = null;
// ── Contribution Bar Chart ──
// ── SHAP Feature Impact Profile Definitions ──
var SHAP_KEYS = [
    { rank: 1, bit: 166, name: "FRAGMENTS", desc: "Number of disconnected molecular fragments", effect: "Decreases" },
    { rank: 2, bit: 99, name: "C=C", desc: "Carbon-carbon double bond", effect: "Decreases" },
    { rank: 3, bit: 44, name: "OTHER", desc: "Contains uncommon/other atom types", effect: "Decreases" },
    { rank: 4, bit: 133, name: "A$A!N", desc: "Ring-chain nitrogen junction", effect: "Increases" },
    { rank: 5, bit: 135, name: "Nnot%A%A", desc: "Non-aromatic nitrogen attached to aromatic system", effect: "Increases" },
    { rank: 6, bit: 35, name: "GROUP IA", desc: "Li, Na, K, Rb, Cs containing structures", effect: "Decreases" },
    { rank: 7, bit: 123, name: "OCO", desc: "Ester/carbonate-like O-C-O motif", effect: "Decreases" },
    { rank: 8, bit: 154, name: "C=O", desc: "Carbonyl group", effect: "Increases" },
    { rank: 9, bit: 139, name: "OH", desc: "Hydroxyl group", effect: "Increases" },
    { rank: 10, bit: 97, name: "NAAAO", desc: "Nitrogen connected through three atoms to oxygen", effect: "Increases" },
    { rank: 11, bit: 65, name: "C%N", desc: "Aromatic carbon-nitrogen bond", effect: "Increases" },
    { rank: 12, bit: 37, name: "NC(O)N", desc: "Urea-like motif", effect: "Increases" },
    { rank: 13, bit: 146, name: "O > 2", desc: "More than two oxygen atoms", effect: "Decreases" },
    { rank: 14, bit: 138, name: "QCH2A > 1", desc: "Multiple heteroatom-CH2-atom motifs", effect: "Increases" },
    { rank: 15, bit: 156, name: "NA(A)A", desc: "Tertiary nitrogen environment", effect: "Increases" },
    { rank: 16, bit: 12, name: "GROUP IB/IIB", desc: "Transition metal-containing structures", effect: "Decreases" },
    { rank: 17, bit: 76, name: "C=C(A)A", desc: "Substituted alkene", effect: "Decreases" },
    { rank: 18, bit: 164, name: "O", desc: "Oxygen atom present", effect: "Increases" },
    { rank: 19, bit: 106, name: "QA(Q)Q", desc: "Atom bonded to multiple heteroatoms", effect: "Increases" },
    { rank: 20, bit: 89, name: "OAAAO", desc: "Oxygen separated by three atoms from oxygen", effect: "Decreases" }
];

function renderFingerprintDonut(chartData, predictionLabel) {
    var mappedData = SHAP_KEYS.map(function (key) {

        var match = chartData.find(function (d) {
            return d.bit_index === key.bit;
        });

        return {
            isPresent: match ? match.status === 'Present' : false,
            direction: FEATURE_DIRECTION[key.bit] || "NON_HERBICIDE"
        };
    });

    var supportsPrediction = [];
    var opposesPrediction = [];

    mappedData.forEach(function (d) {

        var supportsHerbicideSignal =
            (d.direction === "HERBICIDE" && d.isPresent) ||
            (d.direction === "NON_HERBICIDE" && !d.isPresent);

        if (predictionLabel === "Herbicide") {

            if (supportsHerbicideSignal) {
                supportsPrediction.push(d);
            } else {
                opposesPrediction.push(d);
            }

        } else {

            if (!supportsHerbicideSignal) {
                supportsPrediction.push(d);
            } else {
                opposesPrediction.push(d);
            }
        }
    });

    var supportCount = supportsPrediction.length;
    var opposeCount = opposesPrediction.length;

    var ctx =
        document
            .getElementById('fingerprintDonutChart')
            .getContext('2d');

    if (fingerprintDonutChart) {
        fingerprintDonutChart.destroy();
    }

    fingerprintDonutChart = new Chart(ctx, {

        type: 'doughnut',

        data: {

            labels: [
                'Supports Prediction',
                'Opposes Prediction'
            ],

            datasets: [{

                data: [
                    supportCount,
                    opposeCount

                ],

                backgroundColor: [
                    '#2E7D32',
                    '#C62828'

                ],

                borderWidth: 0

            }]
        },

        options: {

            responsive: true,

            cutout: '65%',

            plugins: {

                legend: {
                    position: 'bottom'
                },

                tooltip: {

                    callbacks: {

                        label: function (context) {

                            return (
                                context.label +
                                ': ' +
                                context.raw
                            );

                        }
                    }
                }
            }
        }
    });
}

var FEATURE_DIRECTION = {
    166: "NON_HERBICIDE",
    99: "NON_HERBICIDE",
    44: "NON_HERBICIDE",
    133: "HERBICIDE",
    135: "HERBICIDE",
    35: "NON_HERBICIDE",
    123: "NON_HERBICIDE",
    154: "HERBICIDE",
    139: "HERBICIDE",
    97: "HERBICIDE",
    65: "HERBICIDE",
    37: "HERBICIDE",
    146: "NON_HERBICIDE",
    138: "HERBICIDE",
    156: "HERBICIDE",
    12: "NON_HERBICIDE",
    76: "NON_HERBICIDE",
    164: "HERBICIDE",
    106: "HERBICIDE",
    89: "NON_HERBICIDE"
};

function renderEvidenceDashboard(chartData, predictionLabel) {
    var tbody = document.querySelector('#evidence-table tbody');
    tbody.innerHTML = '';

    var mappedData = SHAP_KEYS.map(function (key) {

        var match = chartData.find(function (d) {
            return d.bit_index === key.bit;
        });

        var isPresent = match ? match.status === 'Present' : false;

        var direction = FEATURE_DIRECTION[key.bit] || "NON_HERBICIDE";

        return {
            rank: key.rank,
            bit: key.bit,
            name: key.name,
            desc: key.desc,
            isPresent: isPresent,
            direction: direction,
            association: (direction === "HERBICIDE")
                ? "Herbicide-associated"
                : "Non-Herbicide-associated"
        };
    });

    mappedData.sort(function (a, b) { return a.rank - b.rank; });

    // Populate Table
    mappedData.forEach(function (d) {
        var tr = document.createElement('tr');
        var statusHtml = d.isPresent
            ? '<span style="color:var(--success);font-weight:700;">Present</span>'
            : '<span style="color:var(--text-secondary);">Absent</span>';
        var assocColor = d.association === 'Herbicide-associated' ? 'var(--success)' : 'var(--danger)';
        var assocHtml = '<span style="color:' + assocColor + ';font-weight:600;">' + d.association + '</span>';

        tr.innerHTML =
            '<td>' + d.rank + '</td>' +
            '<td>Key ' + d.bit + ': ' + d.name + '</td>' +
            '<td>' + d.desc + '</td>' +
            '<td>' + statusHtml + '</td>' +
            '<td>' + assocHtml + '</td>';
        tbody.appendChild(tr);
    });

    function generateCategoryHtml(title, items, color) {
        var html = '<div style="margin-bottom: 1rem;">';
        html += '<h4 style="color:' + color + '; margin-bottom: 0.5rem; font-size: 0.9rem;">' + title + ' (' + items.length + ')</h4>';
        if (items.length === 0) {
            html += '<p style="color:var(--text-secondary); font-size: 0.85rem; margin-top: 0;">None detected</p>';
        } else {
            html += '<div style="display:flex; flex-wrap:wrap; gap:0.4rem;">';
            items.forEach(function (item) {
                var pillBg = color === 'var(--success)' ? 'var(--success-bg)' : (color === 'var(--danger)' ? 'var(--danger-bg)' : '#f5f5f5');
                html += '<span class="feature-pill" style="border:1px solid ' + color + '; color: ' + color + '; background: ' + pillBg + ';">#' + item.rank + ' ' + item.name + '</span>';
            });
            html += '</div>';
        }
        html += '</div>';
        return html;
    }
    var supportsPrediction = [];
    var opposesPrediction = [];

    mappedData.forEach(function (d) {

        var supportsHerbicideSignal =
            (d.direction === "HERBICIDE" && d.isPresent) ||
            (d.direction === "NON_HERBICIDE" && !d.isPresent);

        if (predictionLabel === "Herbicide") {

            if (supportsHerbicideSignal) {
                supportsPrediction.push(d);
            } else {
                opposesPrediction.push(d);
            }

        } else {

            // flipped logic for Non-Herbicide
            if (!supportsHerbicideSignal) {
                supportsPrediction.push(d);
            } else {
                opposesPrediction.push(d);
            }
        }
    });
    var container = document.getElementById('evidence-classification-container');

    if (predictionLabel === "Herbicide") {

        container.innerHTML =
            generateCategoryHtml(
                'Supports Herbicide prediction',
                supportsPrediction,
                'var(--success)'
            ) +
            generateCategoryHtml(
                'Opposes Herbicide prediction',
                opposesPrediction,
                'var(--danger)'
            );

    } else {

        container.innerHTML =
            generateCategoryHtml(
                'Supports Non-Herbicide prediction',
                supportsPrediction,
                'var(--danger)'
            ) +
            generateCategoryHtml(
                'Opposes Non-Herbicide prediction',
                opposesPrediction,
                'var(--success)'
            );
    }
}

// ── Radar Chart ──
function renderRadarChart(properties, predictionLabel) {
    var ctx = document.getElementById('radarChart').getContext('2d');
    var isHerbicide = predictionLabel === 'Herbicide';

    var bgColor = isHerbicide ? 'rgba(27, 94, 32, 0.12)' : 'rgba(183, 28, 28, 0.12)';
    var borderColor = isHerbicide ? 'rgba(27, 94, 32, 0.7)' : 'rgba(183, 28, 28, 0.7)';
    var pointColor = isHerbicide ? '#1B5E20' : '#B71C1C';

    if (radarChartInstance) radarChartInstance.destroy();

    radarChartInstance = new Chart(ctx, {
        type: 'radar',
        data: {
            labels: ['Mol Wt', 'LogP', 'H-Donors', 'H-Acceptors', 'TPSA', 'Rot. Bonds'],
            datasets: [{
                label: 'Compound',
                data: [
                    properties.mw / 100,
                    properties.logp,
                    properties.hbd,
                    properties.hba,
                    properties.tpsa / 10,
                    properties.rotatable_bonds
                ],
                backgroundColor: bgColor,
                borderColor: borderColor,
                pointBackgroundColor: pointColor,
                pointBorderColor: '#FFFFFF',
                pointBorderWidth: 2,
                pointRadius: 4,
                borderWidth: 2,
                fill: true
            }]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            scales: {
                r: {
                    angleLines: { color: 'rgba(0,0,0,0.08)' },
                    grid: { color: 'rgba(0,0,0,0.08)' },
                    pointLabels: {
                        color: '#263238',
                        font: { family: 'Inter', size: 11, weight: '500' }
                    },
                    ticks: { display: false },
                    suggestedMin: 0
                }
            },
            plugins: {
                legend: { display: false },
                tooltip: {
                    backgroundColor: '#263238',
                    cornerRadius: 6,
                    padding: 10,
                    callbacks: {
                        label: function (context) {
                            var val = context.raw;
                            if (context.dataIndex === 0) val = val * 100;
                            if (context.dataIndex === 4) val = val * 10;
                            return ' ' + (Math.round(val * 100) / 100);
                        }
                    }
                }
            }
        }
    });
}

// ── Property Cards ──
function renderPropertyCards(props) {
    var grid = document.getElementById('properties-grid');
    var items = [
        { label: 'Mol. Weight', value: props.mw, unit: 'g/mol' },
        { label: 'LogP', value: props.logp, unit: '' },
        { label: 'H-Donors', value: props.hbd, unit: '' },
        { label: 'H-Acceptors', value: props.hba, unit: '' },
        { label: 'TPSA', value: props.tpsa, unit: 'Å²' },
        { label: 'Rot. Bonds', value: props.rotatable_bonds, unit: '' }
    ];

    grid.innerHTML = items.map(function (i) {
        var unitHtml = i.unit ? '<span style="font-size:0.6em;color:var(--text-secondary);margin-left:2px">' + i.unit + '</span>' : '';
        return '<div class="property-card">' +
            '<div class="property-value">' + i.value + unitHtml + '</div>' +
            '<div class="property-label">' + i.label + '</div>' +
            '</div>';
    }).join('');
}

// ── Single Prediction ──
async function predictSingle() {
    var smiles = document.getElementById('smiles-input').value.trim();
    if (!smiles) {
        showToast('Enter a SMILES string.', 'info');
        return;
    }

    var loading = document.getElementById('loading-single');
    var resultBox = document.getElementById('result-single');

    loading.classList.add('active');
    resultBox.style.display = 'none';

    try {
        var response = await fetch('/herb-pred/predict', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ smiles: smiles })
        });

        var data = await response.json();
        if (!response.ok) throw new Error(data.error || 'Prediction failed');

        // Badge
        var badge = document.getElementById('pred-badge');
        badge.textContent = data.prediction;
        badge.className = 'result-badge ' + (data.prediction === 'Herbicide' ? 'herbicide' : 'not-herbicide');

        // Confidence
        var confValue = parseFloat(data.confidence);
        document.getElementById('pred-conf').textContent = data.confidence;
        var bar = document.getElementById('confidence-bar');
        bar.style.width = '0%';
        requestAnimationFrame(function () { bar.style.width = confValue + '%'; });

        // Properties
        renderPropertyCards(data.properties);

        // Explanation charts
        if (data.explanation && data.explanation.chart_data) {
            renderFingerprintDonut(
                data.explanation.chart_data,
                data.prediction
            );

            renderEvidenceDashboard(
                data.explanation.chart_data,
                data.prediction
            );
        }

        // Radar
        renderRadarChart(data.properties, data.prediction);

        resultBox.style.display = 'block';
        setTimeout(function () {
            resultBox.scrollIntoView({ behavior: 'smooth', block: 'start' });
        }, 100);

    } catch (error) {
        showToast('Error: ' + error.message);
    } finally {
        loading.classList.remove('active');
    }
}

// ── Batch Prediction ──
var currentBatchData = [];

async function predictBatch() {
    var inputText = document.getElementById('batch-input').value;
    var smilesList = inputText.split('\n').map(function (s) { return s.trim(); }).filter(function (s) { return s; });

    if (smilesList.length === 0) {
        showToast('Enter at least one SMILES string.', 'info');
        return;
    }

    var loading = document.getElementById('loading-batch');
    var resultBox = document.getElementById('result-batch');
    var tbody = document.querySelector('#batch-table tbody');

    loading.classList.add('active');
    resultBox.style.display = 'none';
    tbody.innerHTML = '';
    currentBatchData = [];

    try {
        var response = await fetch('/herb-pred/predict_batch', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ smiles_list: smilesList })
        });

        var data = await response.json();
        if (!response.ok) throw new Error(data.error || 'Batch prediction failed');

        currentBatchData = data.results;

        currentBatchData.forEach(function (item) {
            var tr = document.createElement('tr');
            if (item.error) {
                tr.innerHTML =
                    '<td style="font-family:monospace;font-size:0.85rem;">' + item.smiles + '</td>' +
                    '<td colspan="4" style="color:var(--danger)">Error: ' + item.error + '</td>';
            } else {
                var isHerb = item.prediction === 'Herbicide';
                tr.innerHTML =
                    '<td style="font-family:monospace;font-size:0.85rem;">' + item.smiles + '</td>' +
                    '<td><span style="font-weight:600;color:' + (isHerb ? 'var(--success)' : 'var(--danger)') + ';">' + item.prediction + '</span></td>' +
                    '<td>' + item.confidence + '</td>' +
                    '<td>' + item.properties.mw + '</td>' +
                    '<td>' + item.properties.logp + '</td>';
            }
            tbody.appendChild(tr);
        });

        resultBox.style.display = 'block';
    } catch (error) {
        showToast('Error: ' + error.message);
    } finally {
        loading.classList.remove('active');
    }
}

// ── Download CSV ──
function downloadCSV() {
    if (currentBatchData.length === 0) return;

    var csv = "SMILES,Prediction,Confidence,MW,LogP,H-Donors,H-Acceptors,TPSA,Rotatable_Bonds\n";

    currentBatchData.forEach(function (item) {
        if (item.error) {
            csv += '"' + item.smiles + '","Error: ' + item.error + '","","","","","","",""\n';
        } else {
            var p = item.properties;
            csv += '"' + item.smiles + '","' + item.prediction + '","' + item.confidence + '","' + p.mw + '","' + p.logp + '","' + p.hbd + '","' + p.hba + '","' + p.tpsa + '","' + p.rotatable_bonds + '"\n';
        }
    });

    var blob = new Blob([csv], { type: 'text/csv;charset=utf-8;' });
    var url = URL.createObjectURL(blob);
    var link = document.createElement('a');
    link.href = url;
    link.download = 'herbicide_predictions.csv';
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    URL.revokeObjectURL(url);
}

// ── Keyboard shortcuts ──
document.addEventListener('keydown', function (e) {
    if (e.key === 'Enter' && !e.shiftKey) {
        if (document.activeElement && document.activeElement.id === 'smiles-input') {
            e.preventDefault();
            predictSingle();
        }
        if (document.activeElement && document.activeElement.id === 'pubchem-cid') {
            e.preventDefault();
            lookupAndPredict();
        }
    }
});

function downloadEvidenceCSV() {
    const table = document.getElementById("evidence-table");

    if (!table) {
        alert("Evidence table not found.");
        return;
    }

    let csv = [];

    // Headers
    const headers = [];
    table.querySelectorAll("thead th").forEach(th => {
        headers.push(`"${th.innerText.trim()}"`);
    });
    csv.push(headers.join(","));

    // Rows
    table.querySelectorAll("tbody tr").forEach(row => {
        const rowData = [];
        row.querySelectorAll("td").forEach(cell => {
            let text = cell.innerText.replace(/"/g, '""').trim();
            rowData.push(`"${text}"`);
        });
        csv.push(rowData.join(","));
    });

    // Create file
    const csvContent = csv.join("\n");
    const blob = new Blob([csvContent], {
        type: "text/csv;charset=utf-8;"
    });

    const url = URL.createObjectURL(blob);

    const link = document.createElement("a");
    link.href = url;
    link.download = "single_prediction_evidence.csv";

    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);

    URL.revokeObjectURL(url);
}

function insertExampleSmiles() {
    document.getElementById("smiles-input").value =
        "CC1=C(C=CC(=C1Cl)OC(C)C(=O)NC2=CC=CC=C2)Cl";
}