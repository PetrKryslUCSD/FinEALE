classdef  Abaqus_lck_watcher  < handle
    %     Simple class to wait for and watch the lock  (LCK) file.
    %
   
    
    methods % constructor
        
        function self = Abaqus_lck_watcher(Parameters)
            if nargin <1
                Parameters = struct([]);
            end
            
        end
        
    end
    
    
    methods % concrete methods
        
        function A= wait(self,File)
            while true
                if (exist(File))
                    break
                end
            end
            while true
                if (~exist(File))
                    break
                end
            end
        end
        
    end
    
end


