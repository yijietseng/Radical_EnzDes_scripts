import pickle as pkl

def save(arr,fileName):
    fileObject = open(fileName, 'wb')
    pkl.dump(arr, fileObject)
    fileObject.close()
def load(fileName):
    fileObject2 = open(fileName, 'rb')
    modelInput = pkl.load(fileObject2)
    fileObject2.close()
    return modelInput


hbDesign = ["repeat 5",
			"scale:fa_rep 0.1",
			'weight:hbond_sc 4', 
			'weight:hbond_bb_sc 4',
			'weight:fa_intra_rep 0.2',
			"repack",
			"scale:fa_rep 0.1",
			'weight:hbond_sc 4', 
			'weight:hbond_bb_sc 4',
			'weight:fa_intra_rep 0.2',
			"min 0.01",
			"scale:fa_rep 0.265",
			'weight:hbond_sc 4', 
			'weight:hbond_bb_sc 4',
			'weight:fa_intra_rep 0.2',
			"repack",
			"scale:fa_rep 0.280",
			'weight:hbond_sc 4', 
			'weight:hbond_bb_sc 4',
			'weight:fa_intra_rep 0.2',
			"min 0.01",
			"scale:fa_rep 0.559",
			'weight:hbond_sc 4', 
			'weight:hbond_bb_sc 4',
			'weight:fa_intra_rep 0.2',
			"repack",
			"scale:fa_rep 0.581",
			'weight:hbond_sc 4', 
			'weight:hbond_bb_sc 4',
			'weight:fa_intra_rep 0.2',
			"min 0.01",
			"scale:fa_rep 1",
			'weight:hbond_sc 4', 
			'weight:hbond_bb_sc 4',
			'weight:fa_intra_rep 0.2',
			"repack",
			"min 0.00001",
			"accept_to_best",
			"endrepeat"]

hpDesign = ['repeat 5', 
			'scale:fa_rep 0.1', 
			'weight:hbond_sc 4', 
			'weight:hbond_bb_sc 4',
			'weight:fa_intra_rep 0.2',
			'repack',
			'scale:fa_rep 0.1',
			'weight:hbond_sc 4',
			'weight:hbond_bb_sc 4',
			'weight:fa_intra_rep 0.2',
			'min 0.01',
			'scale:fa_rep 0.265',
			'weight:hbond_sc 3',
			'weight:hbond_bb_sc 3',
			'weight:fa_intra_rep 0.15',
			'repack',
			'scale:fa_rep 0.280',
			'weight:hbond_sc 3',
			'weight:hbond_bb_sc 3',
			'weight:fa_intra_rep 0.15',
			'min 0.01',
			'scale:fa_rep 0.559',
			'weight:hbond_sc 2',
			'weight:hbond_bb_sc 2',
			'weight:fa_intra_rep 0.1',
			'repack',
			'scale:fa_rep 0.581',
			'weight:hbond_sc 2',
			'weight:hbond_bb_sc 2',
			'weight:fa_intra_rep 0.1',
			'min 0.01',
			'scale:fa_rep 1',
			'weight:hbond_sc 1',
			'weight:hbond_bb_sc 1',
			'weight:fa_intra_rep 0.005',
			'repack',
			'min 0.00001',
			'accept_to_best',
			'endrepeat']

save(hbDesign, 'hb_RelaxScript_v2.pkl')
save(hpDesign, 'hp_RelaxScript_v2.pkl')
