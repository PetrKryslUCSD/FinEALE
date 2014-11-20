function run_verification
% Run all verification problems.
% 
% function run_verification
% 
% The examples in the folders from the verification_list (see below) are
% all run automatically. The results of the verification suite are recorded
% in a .log file in the working directory.
% 
% The function run_all_examples is invoked to carry out the computation.
% 
    verification_list{1}=[fineale_path filesep 'verification' ];
    run_all_examples(verification_list);
end