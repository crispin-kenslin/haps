from flask import Flask, request, jsonify, render_template
import joblib, os, numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, MACCSkeys, Lipinski
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
from waitress import serve

app = Flask(__name__, static_url_path='/herb-pred/static')

# Load Models
model = joblib.load("models/MACCS_SVM.pkl")
selector = joblib.load("models/MACCS_selector.pkl")
morgan_gen = GetMorganGenerator(radius=2, fpSize=2048)

# Standard 166 public MACCS Key definitions
MACCS_DESCRIPTIONS = {
    1: "Isotope", 2: "103 < Atomic No. < 256", 3: "Group IVA,VA,VIA Periods 4-6", 4: "Actinide",
    5: "Group IIIB,IVB (Sc...)", 6: "Lanthanide", 7: "Group VB,VIB,VIIB", 8: "QAAA@1",
    9: "Group VIII (Fe...)", 10: "Group IIA (Alkaline Earth)", 11: "4-Membered Ring", 12: "Group IB,IIB (Cu...)",
    13: "ON(C)C", 14: "S-S Bond", 15: "OC(O)O", 16: "QAA@1", 17: "CTC", 18: "Group IIIA (B...)",
    19: "7-Membered Ring", 20: "Silicon (Si)", 21: "C=C(Q)Q", 22: "3-Membered Ring", 23: "NC(O)O",
    24: "N-O Bond", 25: "NC(N)N", 26: "C$=C($A)$A", 27: "Iodine (I)", 28: "QCH2Q",
    29: "Phosphorus (P)", 30: "CQ(C)(C)A", 31: "QX", 32: "CSN", 33: "NS", 34: "CH2=A",
    35: "Group IA (Alkali Metal)", 36: "S Heterocycle", 37: "NC(O)N", 38: "NC(C)N", 39: "OS(O)O",
    40: "S-O Bond", 41: "CTN", 42: "Fluorine (F)", 43: "QHAQH", 44: "Other Atoms", 45: "C=CN",
    46: "Bromine (Br)", 47: "SAN", 48: "OQ(O)O", 49: "Charge", 50: "C=C(C)C", 51: "CSO",
    52: "NN Bond", 53: "QHAAAQH", 54: "QHAAQH", 55: "OSO", 56: "ON(O)C", 57: "O Heterocycle",
    58: "QSQ", 59: "Snot%A%A", 60: "S=O Bond", 61: "AS(A)A", 62: "A$A!A$A", 63: "N=O Bond",
    64: "A$A!S", 65: "C%N (Nitrile)", 66: "CC(C)(C)A", 67: "QS", 68: "QHQH", 69: "QQH",
    70: "QNQ", 71: "NO Bond", 72: "OAAO", 73: "S=A", 74: "CH3ACH3", 75: "A!N$A",
    76: "C=C(A)A", 77: "NAN", 78: "C=N Bond", 79: "NAAN", 80: "NAAAN", 81: "SA(A)A",
    82: "ACH2QH", 83: "QAAAA@1", 84: "NH2 Group", 85: "CN(C)C", 86: "CH2QCH2", 87: "X!A$A",
    88: "Sulfur Atom (S)", 89: "OAAAO", 90: "QHAACH2A", 91: "QHAAACH2A", 92: "OC(N)C",
    93: "QCH3", 94: "QN", 95: "NAAO", 96: "5-Membered Ring", 97: "N NAAAO / NH2", 98: "QAAAAA@1",
    99: "C=C Double Bond", 100: "ACH2N", 101: "8-Membered Ring", 102: "QO", 103: "Chlorine (Cl)",
    104: "QHACH2A", 105: "A$A($A)$A", 106: "QA(Q)Q", 107: "XA(A)A", 108: "CH3AAACH2A",
    109: "ACH2O", 110: "NCO", 111: "NACH2A", 112: "AA(A)(A)A", 113: "Onot%A%A",
    114: "CH3CH2A", 115: "CH3ACH2A", 116: "CH3AACH2A", 117: "NAO", 118: "ACH2CH2A > 1",
    119: "N=A", 120: "Heterocyclic Atom > 1", 121: "N Heterocycle", 122: "AN(A)A", 123: "OCO",
    124: "QQ", 125: "Aromatic Ring > 1", 126: "A!O!A", 127: "A$A!O > 1", 128: "ACH2AAACH2A",
    129: "ACH2AACH2A", 130: "QQ > 1", 131: "QH > 1", 132: "OACH2A", 133: "A$A!N",
    134: "Halogen Atom (X)", 135: "Nnot%A%A", 136: "O=A > 1", 137: "Heterocycle present",
    138: "QCH2A > 1", 139: "OH (Hydroxyl Group)", 140: "Oxygen Atoms > 3", 141: "Methyl Groups > 2",
    142: "Nitrogen Atoms > 1", 143: "A$A!O", 144: "Anot%A%Anot%A", 145: "6-Membered Rings > 1",
    146: "Oxygen Atoms > 2", 147: "ACH2CH2A", 148: "AQ(A)A", 149: "Methyl Groups > 1",
    150: "A!A$A!A", 151: "NH Group", 152: "OC(C)C", 153: "QCH2A", 154: "C=O (Carbonyl Group)",
    155: "A!CH2!A", 156: "NA(A)A", 157: "C-O Single Bond", 158: "C-N Single Bond",
    159: "Oxygen Atoms > 1", 160: "CH3 (Methyl Group)", 161: "Nitrogen Atom (N)",
    162: "Aromatic Ring present", 163: "6-Membered Ring present", 164: "Oxygen Atom (O)",
    165: "Ring present", 166: "Fragments present"
}

def get_properties(mol):
    return {
        "mw": float(f"{Descriptors.MolWt(mol):.2f}"),
        "logp": float(f"{Descriptors.MolLogP(mol):.2f}"),
        "hbd": Lipinski.NumHDonors(mol),
        "hba": Lipinski.NumHAcceptors(mol),
        "tpsa": float(f"{Descriptors.TPSA(mol):.2f}"),
        "rotatable_bonds": Lipinski.NumRotatableBonds(mol)
    }

def get_prediction(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol: return None, {"error": "Invalid SMILES"}
        
        # MACCS vector yields 167 bits (index 0 is padding)
        fp = np.array(MACCSkeys.GenMACCSKeys(mol))
        
        # Unmask original MACCS key indexes via the Selector model
        if hasattr(selector, 'get_support'):
            selected_indices = selector.get_support(indices=True)
        else:
            selected_indices = list(range(167))
            
        fp_selected = selector.transform([fp])
        
        pred = model.predict(fp_selected)[0]
        conf = float(np.max(model.predict_proba(fp_selected))) * 100
        
        # Extract feature contributions (supports both Linear and non-linear fallback)
        contributions = {}
        if hasattr(model, 'coef_'):
            # Linear SVM gives precise weights for each selected key
            weights = model.coef_[0]
            for idx, orig_idx in enumerate(selected_indices):
                is_present = int(fp[orig_idx])
                contributions[int(orig_idx)] = float(weights[idx] * is_present)
        elif hasattr(selector, 'scores_'):
            # Fallback to feature selection ANOVA/Chi2 scores if SVM is RBF kernel
            for idx, orig_idx in enumerate(selected_indices):
                if int(fp[orig_idx]) == 1:
                    contributions[int(orig_idx)] = float(selector.scores_[orig_idx])
                else:
                    contributions[int(orig_idx)] = 0.0
        else:
            # Fallback to binary active status
            for orig_idx in selected_indices:
                contributions[int(orig_idx)] = float(fp[orig_idx])

        # Formulate structured chart data for the UI
        chart_data = []
        for orig_idx in selected_indices:
            orig_idx = int(orig_idx)
            if orig_idx == 0: 
                continue # skip filler bit
                
            name = MACCS_DESCRIPTIONS.get(orig_idx, f"Key {orig_idx}")
            is_present = int(fp[orig_idx])
            
            chart_data.append({
                "bit_index": orig_idx,
                "name": name,
                "status": "Present" if is_present else "Absent",
                "value": is_present,
                "contribution": contributions.get(orig_idx, 0.0)
            })
            
        # Sort features so the most impactful fragments appear first in charts
        chart_data = sorted(chart_data, key=lambda x: abs(x["contribution"]), reverse=True)
        
        props = get_properties(mol)
        pred_label = "Herbicide" if pred == 1 else "Not Herbicide"
        
        return {
            "smiles": smiles,
            "prediction": pred_label,
            "confidence": f"{conf:.2f}%",
            "properties": props,
            "explanation": {
                "chart_data": chart_data # Directly useable by front-end graphing engines
            }
        }, None
    except Exception as e:
        import traceback
        print(traceback.format_exc())
        return None, {"error": f"Prediction error: {str(e)}"}

@app.route('/herb-pred/')
def index():
    return render_template('index.html')
MACCS_KEYS = [
(1, "ISOTOPE", "Atom with isotopic labeling"),
(2, "103 < ATOMIC NO < 256", "Very heavy element present"),
(3, "GROUP IVA/VA/VIA (P4-6)", "Heavy p-block element"),
(4, "ACTINIDE", "Actinide metal atom"),
(5, "GROUP IIIB/IVB", "Early transition metal"),
(6, "LANTHANIDE", "Rare earth element"),
(7, "GROUP VB/VIB/VIIB", "Mid transition metal"),
(8, "QAAA@1", "Specific 3-connected ring atom"),
(9, "GROUP VIII", "Fe/Co/Ni-type metals"),
(10, "GROUP IIA", "Alkaline earth metals"),

(11, "4M RING", "Four-membered ring"),
(12, "GROUP IB/IIB", "Cu/Zn metals"),
(13, "ON(C)C", "N bonded to O and C"),
(14, "S-S", "Disulfide bond"),
(15, "OC(O)O", "Carbonate/ester group"),
(16, "QAA@1", "Ring-connected atom"),
(17, "CTC", "Carbon chain linkage"),
(18, "GROUP IIIA", "Boron group element"),
(19, "7M RING", "Seven-membered ring"),
(20, "SI", "Silicon atom"),

(21, "C=C(Q)Q", "Substituted double bond"),
(22, "3M RING", "Three-membered ring"),
(23, "NC(O)O", "Carbamate/urea group"),
(24, "N-O", "Nitrogen-oxygen bond"),
(25, "NC(N)N", "Guanidine-like group"),
(26, "C$=C($A)$A", "Conjugated double bond system"),
(27, "I", "Iodine atom"),
(28, "QCH2Q", "Heteroatom-CH2-heteroatom"),
(29, "P", "Phosphorus atom"),
(30, "CQ(C)(C)A", "Branched carbon center"),

(31, "QX", "Halogen/heteroatom pattern"),
(32, "CSN", "C-S-N fragment"),
(33, "NS", "N-S bond"),
(34, "CH2=A", "Methylene double bond"),
(35, "ALKALI METAL", "Group IA metal"),
(36, "S HETEROCYCLE", "Sulfur-containing ring"),
(37, "NC(O)N", "Urea-like structure"),
(38, "NC(C)N", "Amine-rich structure"),
(39, "OS(O)O", "Oxidized sulfur group"),
(40, "S-O", "Sulfur-oxygen bond"),

(41, "CTN", "Carbon-nitrogen linkage"),
(42, "F", "Fluorine atom"),
(43, "QHAQH", "Hydrogen bonding motif"),
(44, "OTHER", "Unclassified fragment"),
(45, "C=CN", "Alkene attached to nitrogen"),
(46, "BR", "Bromine atom"),
(47, "SAN", "Sulfur-amine network"),
(48, "OQ(O)O", "Oxygen-rich group"),
(49, "CHARGE", "Charged structure"),
(50, "C=C(C)C", "Substituted alkene"),

(51, "CSO", "Sulfoxide-like group"),
(52, "NN", "N-N bond"),
(53, "QHAAAQH", "Hydrogen donor pattern"),
(54, "QHAAQH", "Hydrogen bonding system"),
(55, "OSO", "Sulfoxide/sulfone motif"),
(56, "ON(O)C", "Nitro/amine oxygen system"),
(57, "O HETEROCYCLE", "Oxygen in ring"),
(58, "QSQ", "Sulfur-centered pattern"),
(59, "Snot%A%A", "Complex sulfur fragment"),
(60, "S=O", "Sulfoxide group"),

(61, "AS(A)A", "Aromatic sulfur system"),
(62, "A$A!A$A", "Aromatic ring connectivity"),
(63, "N=O", "Nitroso group"),
(64, "A$A!S", "Aromatic-sulfur linkage"),
(65, "C%N", "Cyanide-like group"),
(66, "CC(C)(C)A", "Branched alkyl carbon"),
(67, "QS", "Sulfur environment"),
(68, "QHQH", "Hydrogen bond pair"),
(69, "QQH", "Heteroatom hydrogen motif"),
(70, "QNQ", "Nitrogen bridge"),

(71, "NO", "Nitrogen-oxygen bond"),
(72, "OAAO", "Oxygen separated by carbon"),
(73, "S=A", "Double-bonded sulfur"),
(74, "CH3ACH3", "Two methyl groups"),
(75, "A!N$A", "Aromatic nitrogen linkage"),
(76, "C=C(A)A", "Substituted alkene"),
(77, "NAN", "Nitrogen chain"),
(78, "C=N", "Imine group"),
(79, "NAAN", "Nitrogen-rich chain"),
(80, "NAAAN", "Extended nitrogen chain"),

(81, "SA(A)A", "Aromatic sulfur system"),
(82, "ACH2QH", "H-bonding methylene"),
(83, "QAAAA@1", "Ring-connected atom"),
(84, "NH2", "Amine group"),
(85, "CN(C)C", "Tertiary amine"),
(86, "CH2QCH2", "Heteroatom-CH2-heteroatom"),
(87, "X!A$A", "Halogen-aromatic pattern"),
(88, "S", "Sulfur atom"),
(89, "OAAAO", "Oxygen chain"),
(90, "QHAACH2A", "H-bond donor chain"),
(91, "QHAAACH2A", "Hydrogen bond donor chain with carbon spacer"),
(92, "OC(N)C", "Oxygen-nitrogen-carbon functional group"),
(93, "QCH3", "Methyl group attached to heteroatom"),
(94, "QN", "Heteroatom linked to nitrogen"),
(95, "NAAO", "Nitrogen-oxygen separated by carbon"),
(96, "5M RING", "Five-membered ring system"),
(97, "NAAAO", "Nitrogen chain ending in oxygen"),
(98, "QAAAAA@1", "Extended ring-connected atom system"),
(99, "C=C", "Carbon-carbon double bond"),
(100, "ACH2N", "Carbon chain linked to nitrogen"),

(101, "8M RING", "Eight-membered ring system"),
(102, "QO", "Heteroatom oxygen environment"),
(103, "CL", "Chlorine atom present"),
(104, "QHACH2A", "Hydrogen bonding through methylene"),
(105, "A$A($A)$A", "Complex aromatic substitution pattern"),
(106, "QA(Q)Q", "Heteroatom connected to multiple groups"),
(107, "XA(A)A", "Halogen attached to aromatic system"),
(108, "CH3AAACH2A", "Methyl group connected via carbon chain"),
(109, "ACH2O", "Carbon chain linked to oxygen"),
(110, "NCO", "Nitrogen-carbon-oxygen functional group"),

(111, "NACH2A", "Nitrogen connected via methylene"),
(112, "AA(A)(A)A", "Highly substituted carbon center"),
(113, "Onot%A%A", "Oxygen in complex aromatic environment"),
(114, "CH3CH2A", "Ethyl group in structure"),
(115, "CH3ACH2A", "Methyl-methylene fragment"),
(116, "CH3AACH2A", "Branched alkyl chain"),
(117, "NAO", "Nitrogen-oxygen linkage"),
(118, "ACH2CH2A", "Carbon chain with two methylenes"),
(119, "N=A", "Imine or double-bonded nitrogen"),
(120, "HETEROCYCLIC ATOM > 1", "Multiple heteroatoms in ring"),

(121, "N HETEROCYCLE", "Nitrogen-containing ring"),
(122, "AN(A)A", "Substituted amine environment"),
(123, "OCO", "Carbonyl-oxygen-carbon group"),
(124, "QQ", "Heteroatom-rich region"),
(125, "AROMATIC RING > 1", "Multiple aromatic rings"),
(126, "A!O!A", "Oxygen linking two carbon atoms"),
(127, "A$A!O > 1", "Multiple aromatic oxygen linkages"),
(128, "ACH2AAACH2A", "Extended carbon chain"),
(129, "ACH2AACH2A", "Branched alkyl chain"),
(130, "QQ > 1", "Multiple heteroatom sites"),

(131, "QH > 1", "Multiple hydrogen bond donors"),
(132, "OACH2A", "Oxygen connected via carbon"),
(133, "A$A!N", "Aromatic nitrogen linkage"),
(134, "X (HALOGEN)", "Halogen atom present"),
(135, "Nnot%A%A", "Nitrogen in complex environment"),
(136, "O=A > 1", "Multiple carbonyl groups"),
(137, "HETEROCYCLE", "Ring with heteroatom"),
(138, "QCH2A > 1", "Multiple heteroatom-linked carbons"),
(139, "OH", "Hydroxyl group"),
(140, "O > 3", "Multiple oxygen atoms"),

(141, "CH3 > 2", "Multiple methyl groups"),
(142, "N > 1", "Multiple nitrogen atoms"),
(143, "A$A!O", "Aromatic oxygen linkage"),
(144, "Anot%A%Anot%A", "Complex aromatic system"),
(145, "6M RING > 1", "Multiple six-membered rings"),
(146, "O > 2", "At least two oxygen atoms"),
(147, "ACH2CH2A", "Ethylene chain fragment"),
(148, "AQ(A)A", "Substituted aromatic center"),
(149, "CH3 > 1", "At least one methyl group"),
(150, "A!A$A!A", "Aromatic connectivity pattern"),
(151, "NH", "Amine hydrogen"),
(152, "OC(C)C", "Alcohol/ether group"),
(153, "QCH2A", "Heteroatom-linked carbon"),
(154, "C=O", "Carbonyl group"),

(155, "A!CH2!A", "Aromatic-CH2-aromatic linkage"),
(156, "NA(A)A", "Nitrogen substitution"),
(157, "C-O", "Alcohol/ether bond"),
(158, "C-N", "Amine bond"),
(159, "O > 1", "Multiple oxygen atoms"),
(160, "CH3", "Methyl group"),

(161, "N", "Nitrogen atom"),
(162, "AROMATIC", "Aromatic system"),
(163, "6M RING", "Six-membered ring"),
(164, "O", "Oxygen atom"),
(165, "RING", "Any cyclic structure"),
(166, "FRAGMENTS", "General molecular fragments"),
]
@app.route('/herb-pred/keys')
def keys():
    return render_template('keys.html', keys = MACCS_KEYS)

@app.route('/herb-pred/about')
def about():
    return render_template('about.html')

@app.route('/herb-pred/help')
def help_page():
    return render_template('help.html')

@app.route('/herb-pred/pubchem/<int:cid>')
def pubchem_lookup(cid):
    """Proxy endpoint to fetch SMILES from PubChem by CID."""
    import urllib.request, json as json_mod
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES,IUPACName,MolecularFormula/JSON"
        req = urllib.request.Request(url, headers={"User-Agent": "HerbicideClassifier/1.0"})
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json_mod.loads(resp.read().decode())
        props = data["PropertyTable"]["Properties"][0]
        return jsonify({
            "cid": cid,
            "smiles": props.get("ConnectivitySMILES", ""),
            "name": props.get("IUPACName", "Unknown"),
            "formula": props.get("MolecularFormula", "")
        })
    except Exception as e:
        return jsonify({"error": f"PubChem lookup failed: {str(e)}"}), 400

@app.route('/herb-pred/predict', methods=['POST'])
def predict():
    smiles = request.json.get("smiles")
    if not smiles: return jsonify({"error": "No SMILES provided"}), 400
    res, err = get_prediction(smiles)
    if err: return jsonify(err), 400
    return jsonify(res)

@app.route('/herb-pred/predict_batch', methods=['POST'])
def predict_batch():
    smiles_list = request.json.get("smiles_list", [])
    if not smiles_list: return jsonify({"error": "No SMILES list provided"}), 400
    
    results = []
    for s in smiles_list:
        s = s.strip()
        if not s: continue
        res, err = get_prediction(s)
        if err:
            results.append({"smiles": s, "error": err["error"]})
        else:
            results.append(res)
            
    return jsonify({"results": results})

if __name__ == '__main__':
    app.run(host="0.0.0.0", port=5000, debug=True, use_reloader=True) 