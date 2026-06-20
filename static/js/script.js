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
var contributionChartInstance = null;
var presenceChartInstance = null;
var radarChartInstance = null;

// ── Contribution Bar Chart ──
function renderContributionChart(chartData) {
    var ctx = document.getElementById('contributionChart').getContext('2d');

    var top = chartData
        .filter(function (d) { return d.contribution !== 0; })
        .sort(function (a, b) { return Math.abs(b.contribution) - Math.abs(a.contribution); })
        .slice(0, 15);

    var labels = top.map(function (d) { return d.name; });
    var values = top.map(function (d) { return +d.contribution.toFixed(3); });
    var colors = values.map(function (v) {
        return v > 0 ? 'rgba(46, 125, 50, 0.75)' : 'rgba(198, 40, 40, 0.75)';
    });

    if (contributionChartInstance) contributionChartInstance.destroy();

    contributionChartInstance = new Chart(ctx, {
        type: 'bar',
        data: {
            labels: labels,
            datasets: [{
                label: 'Contribution',
                data: values,
                backgroundColor: colors,
                borderRadius: 3,
                barPercentage: 0.7,
                categoryPercentage: 0.85
            }]
        },
        options: {
            indexAxis: 'y',
            responsive: true,
            maintainAspectRatio: false,
            layout: { padding: { left: 4, right: 12 } },
            scales: {
                x: {
                    grid: { color: 'rgba(0,0,0,0.06)' },
                    ticks: { color: '#546E7A', font: { family: 'Inter', size: 11 } },
                    title: {
                        display: true,
                        text: 'SVM Weight × Presence',
                        color: '#78909C',
                        font: { family: 'Inter', size: 11 }
                    }
                },
                y: {
                    grid: { display: false },
                    ticks: { color: '#263238', font: { family: 'Inter', size: 15 } }
                }
            },
            plugins: {
                legend: { display: false },
                tooltip: {
                    backgroundColor: '#263238',
                    cornerRadius: 6,
                    padding: 10,
                    titleFont: { family: 'Inter', weight: '600' },
                    bodyFont: { family: 'Inter' },
                    callbacks: {
                        title: function (items) {
                            var idx = items[0].dataIndex;
                            return top[idx].name + ' (Key ' + top[idx].bit_index + ')';
                        },
                        label: function (item) {
                            var d = top[item.dataIndex];
                            return ['Contribution: ' + d.contribution.toFixed(4), 'Status: ' + d.status];
                        }
                    }
                }
            }
        }
    });
}

// ── Presence Donut Chart ──
function renderPresenceChart(chartData) {
    var ctx = document.getElementById('presenceChart').getContext('2d');

    var present = chartData.filter(function (d) { return d.status === 'Present'; }).length;
    var absent = chartData.filter(function (d) { return d.status === 'Absent'; }).length;

    if (presenceChartInstance) presenceChartInstance.destroy();

    presenceChartInstance = new Chart(ctx, {
        type: 'doughnut',
        data: {
            labels: ['Present', 'Absent'],
            datasets: [{
                data: [present, absent],
                backgroundColor: ['rgba(46, 125, 50, 0.7)', 'rgba(207, 216, 220, 0.7)'],
                borderColor: ['#2E7D32', '#B0BEC5'],
                borderWidth: 1.5,
                hoverOffset: 6
            }]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            cutout: '65%',
            plugins: {
                legend: {
                    position: 'bottom',
                    labels: {
                        color: '#263238',
                        font: { family: 'Inter', size: 12 },
                        padding: 16,
                        usePointStyle: true,
                        pointStyleWidth: 10
                    }
                },
                tooltip: {
                    backgroundColor: '#263238',
                    cornerRadius: 6,
                    padding: 10,
                    callbacks: {
                        label: function (item) {
                            var total = present + absent;
                            var pct = ((item.raw / total) * 100).toFixed(1);
                            return ' ' + item.label + ': ' + item.raw + ' keys (' + pct + '%)';
                        }
                    }
                }
            }
        },
        plugins: [{
            id: 'centerText',
            afterDraw: function (chart) {
                var c = chart.ctx;
                var w = chart.width;
                var top = chart.chartArea.top;
                var bottom = chart.chartArea.bottom;
                var cy = (top + bottom) / 2;
                c.save();
                c.textAlign = 'center';
                c.textBaseline = 'middle';
                c.fillStyle = '#263238';
                c.font = "700 1.5rem 'Inter'";
                c.fillText(present, w / 2, cy - 8);
                c.fillStyle = '#78909C';
                c.font = "500 0.65rem 'Inter'";
                c.fillText('KEYS ACTIVE', w / 2, cy + 12);
                c.restore();
            }
        }]
    });
}

// ── Feature Pills ──
function renderFeaturePills(chartData) {
    var container = document.getElementById('feature-pills');
    container.innerHTML = '';

    var topPresent = chartData
        .filter(function (d) { return d.status === 'Present'; })
        .sort(function (a, b) { return Math.abs(b.contribution) - Math.abs(a.contribution); })
        .slice(0, 10);

    topPresent.forEach(function (d) {
        var pill = document.createElement('span');
        pill.className = 'feature-pill present';
        pill.innerHTML = '<span class="dot"></span>' + d.name;
        pill.title = 'Key ' + d.bit_index + ' | Contribution: ' + d.contribution.toFixed(4);
        container.appendChild(pill);
    });
}

// ── Radar Chart ──
function renderRadarChart(properties) {
    var ctx = document.getElementById('radarChart').getContext('2d');

    if (radarChartInstance) radarChartInstance.destroy();

    radarChartInstance = new Chart(ctx, {
        type: 'radar',
        data: {
            labels: ['Mol Wt (/100)', 'LogP', 'H-Donors', 'H-Acceptors', 'TPSA (/10)', 'Rot. Bonds'],
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
                backgroundColor: 'rgba(27, 94, 32, 0.12)',
                borderColor: 'rgba(27, 94, 32, 0.7)',
                pointBackgroundColor: '#1B5E20',
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
            renderContributionChart(data.explanation.chart_data);
            renderPresenceChart(data.explanation.chart_data);
            renderFeaturePills(data.explanation.chart_data);
        }

        // Radar
        renderRadarChart(data.properties);

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