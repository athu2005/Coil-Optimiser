% Get all blocks in the model
blocks = find_system('sim1');  % replace with actual model name

% List subsystems only
subsystems = find_system('sim1', 'BlockType', 'SubSystem');
