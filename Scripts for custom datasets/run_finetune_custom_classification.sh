python finetune_custom_regression.py \
        --device cuda \
        --batch_size 32  \
        --n_head 12 \
        --n_layer 12 \
        --n_embd 768 \
        --d_dropout 0.1 \
        --dropout 0.1 \
        --lr_start 3e-5 \
        --num_workers 8\
        --max_epochs 500 \
        --num_feats 32 \
        --data_root '../data/ADD NAME OF DATASET FOLDER HERE (MAKE SURE IT'S IN THE DATA FOLDER)" \
        --seed_path '../data/checkpoints/Checkpoint_3_30000.ckpt' \
        --dataset_name "THIS SHOULD BE THE NAME OF YOUR DATASET FOLDER AGAIN" \
        --measure_name "THIS SHOULD BE THE COLUMN NAME OF YOUR TARGET VARIABLE" \
        --dims 768 768 768 1 \
        --checkpoints_folder 'CHOOSE A NAME FOR YOUR CHECKPOINT FOLDER'\
        --num_classes 2 \

